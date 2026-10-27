"""
ANP pipeline inventory: harmonisation, onshore/offshore classification and the
line geometry used for every distance computation.

v3 kept every segment intersecting NE Brazil buffered by 1 degree (~110 km),
which pulled in 28 Espirito Santo segments and 127 mostly-offshore segments,
and reported all 417 as "within the study domain". Here each segment is
classified by where its length actually lies, and distances are measured only
to the ONSHORE parts of the network (retrievals over water are unreliable and
pixels over the sea are excluded anyway).
"""
from pathlib import Path

import geopandas as gpd
import pandas as pd
from pyproj import Geod
from shapely.ops import transform

from .config import Config

# Lengths are geodesic on the WGS84 ellipsoid (a projected CRS such as the
# Brazil Polyconic overstates N-S lengths by 2-5 % this far from its meridian).
GEOD = Geod(ellps="WGS84")


def geodesic_km(geoms):
    return geoms.map(lambda g: GEOD.geometry_length(g) / 1000.0 if g is not None and not g.is_empty else 0.0)


def _force_2d(geom):
    if geom is None:
        return None
    return transform(lambda x, y, z=None: (x, y), geom)


def load_infrastructure(cfg: Config, states: gpd.GeoDataFrame, ne: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """All ANP segments that touch NE land (or lie within 1 degree of it), with
    length shares on NE land, on other Brazilian land and offshore."""
    ne_union = ne.union_all()
    brazil_union = states.union_all()
    window = ne_union.buffer(1.0)

    parts = []
    for fname, (cls, name_col) in cfg.infra_files.items():
        path = Path(cfg.base_dir) / cfg.infra_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"ANP layer not found: {path}")
        g = gpd.read_file(path).to_crs(4326)
        g["asset_name"] = g[name_col].fillna("Not identified").astype(str).str.strip()
        g["asset_class"] = cls
        g["source_file"] = fname
        g = g[["asset_name", "asset_class", "source_file", "geometry"]]
        g = g[g.geometry.notna() & ~g.geometry.is_empty]
        g = g[g.intersects(window)].copy()
        g["geometry"] = g.geometry.map(_force_2d)
        parts.append(g)
    infra = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=4326)
    infra["segment_id"] = [f"S{i:04d}" for i in range(len(infra))]

    total = geodesic_km(infra.geometry)
    on_ne = geodesic_km(infra.geometry.intersection(ne_union))
    on_land = geodesic_km(infra.geometry.intersection(brazil_union))
    infra["length_km"] = total
    infra["length_ne_land_km"] = on_ne
    infra["share_ne_land"] = (on_ne / total).fillna(0)
    infra["share_offshore"] = (1 - on_land / total).clip(lower=0).fillna(0)
    infra["share_other_state"] = ((on_land - on_ne) / total).clip(lower=0).fillna(0)

    def _location(r):
        if r.share_ne_land >= cfg.onshore_share_min:
            return "onshore NE"
        if r.share_offshore >= r.share_other_state:
            return "offshore"
        return "other state"

    infra["location"] = infra.apply(_location, axis=1)
    return infra


def onshore_lines(infra: gpd.GeoDataFrame, states: gpd.GeoDataFrame, cfg: Config) -> gpd.GeoDataFrame:
    """Parts of every segment that lie on Brazilian land (any state), simplified.
    Pixels near the NE border therefore also see pipelines just across it."""
    land = states.union_all()
    lines = infra[["segment_id", "asset_name", "asset_class", "location", "geometry"]].copy()
    lines["geometry"] = lines.geometry.intersection(land)
    lines = lines[~lines.geometry.is_empty & lines.geometry.notna()]
    lines = lines.explode(index_parts=False)
    lines = lines[lines.geom_type == "LineString"]
    lines["geometry"] = lines.geometry.simplify(cfg.simplify_tolerance_deg)
    return lines.reset_index(drop=True)


def inventory_table(infra: gpd.GeoDataFrame) -> pd.DataFrame:
    """Manuscript Table 1: segments and length by class and location."""
    counts = pd.crosstab(infra["asset_class"], infra["location"])
    for col in ("onshore NE", "offshore", "other state"):
        if col not in counts:
            counts[col] = 0
    counts = counts[["onshore NE", "offshore", "other state"]]
    counts.columns = ["segments_onshore_NE", "segments_offshore", "segments_other_state"]
    length = infra.groupby("asset_class")["length_ne_land_km"].sum().rename("length_on_NE_land_km")
    tab = counts.join(length)
    tab.loc["Total"] = tab.sum(numeric_only=True)
    tab["segments_total"] = tab[["segments_onshore_NE", "segments_offshore", "segments_other_state"]].sum(axis=1)
    tab["length_on_NE_land_km"] = tab["length_on_NE_land_km"].round(1)
    cnt = ["segments_onshore_NE", "segments_offshore", "segments_other_state", "segments_total"]
    tab[cnt] = tab[cnt].astype(int)
    return tab.reset_index().rename(columns={"asset_class": "infrastructure_type", "index": "infrastructure_type"})
