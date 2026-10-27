"""
State boundaries and coastline.

State identity is ALWAYS carried as an attribute (UF acronym + name) of the
geometry. The v3 notebook dropped the attributes before sending the polygons
to Earth Engine and then rebuilt the names from a hard-coded index map in IBGE
order, while the source GeoJSON is in alphabetical order: 6 of 9 states got
another state's statistics. Nothing here depends on feature order.
"""
from pathlib import Path

import geopandas as gpd

from .config import Config

_UF_COLUMNS = ("sigla", "SIGLA_UF", "SIGLA", "UF", "abbrev_state")
_NAME_COLUMNS = ("name", "NM_UF", "NM_ESTADO", "nome", "name_state")


def _pick(columns, candidates, what):
    for c in candidates:
        if c in columns:
            return c
    raise KeyError(f"No {what} column found among {list(columns)}")


def load_states(cfg: Config) -> gpd.GeoDataFrame:
    """All Brazilian states with standardised columns ['UF', 'state_name', 'geometry'] in EPSG:4326."""
    src = cfg.states_url
    if cfg.states_file:
        p = Path(cfg.states_file)
        src = p if p.is_absolute() else Path(cfg.base_dir) / p
    gdf = gpd.read_file(src)
    uf_col = _pick(gdf.columns, _UF_COLUMNS, "UF acronym")
    nm_col = _pick(gdf.columns, _NAME_COLUMNS, "state name")
    gdf = gdf.rename(columns={uf_col: "UF", nm_col: "state_name"})[["UF", "state_name", "geometry"]]
    gdf = gdf.to_crs(4326)
    gdf["geometry"] = gdf.geometry.make_valid()
    return gdf


def ne_states(states: gpd.GeoDataFrame, cfg: Config) -> gpd.GeoDataFrame:
    """NE states sorted by UF, with a 1-based integer id used by the state raster."""
    ne = states[states["UF"].isin(cfg.ne_states)].sort_values("UF").reset_index(drop=True)
    missing = set(cfg.ne_states) - set(ne["UF"])
    if missing:
        raise ValueError(f"States missing from boundary file: {sorted(missing)}")
    ne["state_id"] = range(1, len(ne) + 1)
    return ne


def coastline(states: gpd.GeoDataFrame):
    """Boundary of the Brazilian land mass. Within NE Brazil this is the Atlantic
    coastline (international borders are > 1000 km away)."""
    return states.union_all().boundary
