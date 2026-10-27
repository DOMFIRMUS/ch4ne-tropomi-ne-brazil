"""
Raster I/O on the analysis grid and geometry-derived layers (state ids, land
fraction, distances to pipelines and to the coast).

All statistics are computed on the full pixel census of the exported rasters,
so every number in the paper can be recomputed locally and exactly.
"""
import json
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from affine import Affine
from rasterio.features import rasterize
from shapely import STRtree, points

from .config import Config, NODATA, STATS_BANDS

# WGS84 ellipsoid
WGS84_A_KM = 6378.137
WGS84_E2 = 0.00669437999014


def _authalic_q(phi):
    e = np.sqrt(WGS84_E2)
    s = np.sin(phi)
    return (1 - WGS84_E2) * (s / (1 - WGS84_E2 * s ** 2) - np.log((1 - e * s) / (1 + e * s)) / (2 * e))


@dataclass
class Grid:
    transform: Affine
    width: int
    height: int

    @property
    def res(self) -> float:
        return self.transform.a

    @property
    def lon(self) -> np.ndarray:
        """Pixel-centre longitudes, shape (height, width)."""
        cols = np.arange(self.width) + 0.5
        return np.broadcast_to(self.transform.c + cols * self.transform.a, (self.height, self.width))

    @property
    def lat(self) -> np.ndarray:
        rows = np.arange(self.height) + 0.5
        return np.broadcast_to((self.transform.f + rows * self.transform.e)[:, None], (self.height, self.width))

    @property
    def extent(self):
        """matplotlib imshow extent (left, right, bottom, top)."""
        left, top = self.transform.c, self.transform.f
        return (left, left + self.width * self.res, top - self.height * self.res, top)

    @property
    def pixel_area_km2(self) -> np.ndarray:
        """Exact area of each lat/lon cell on the WGS84 ellipsoid."""
        d = np.deg2rad(self.res)
        lat_top = np.deg2rad(self.transform.f + np.arange(self.height) * self.transform.e)
        lat_bot = lat_top + np.deg2rad(self.transform.e)
        a = 0.5 * WGS84_A_KM ** 2 * d * np.abs(_authalic_q(lat_top) - _authalic_q(lat_bot))
        return np.broadcast_to(a[:, None], (self.height, self.width))

    def rowcol(self, lon, lat):
        col = np.floor((np.asarray(lon) - self.transform.c) / self.res).astype(int)
        row = np.floor((self.transform.f - np.asarray(lat)) / self.res).astype(int)
        return row, col


def grid_from_config(cfg: Config) -> Grid:
    w = int(round((cfg.grid_east - cfg.grid_west) / cfg.grid_res_deg))
    h = int(round((cfg.grid_north - cfg.grid_south) / cfg.grid_res_deg))
    t = Affine(cfg.grid_res_deg, 0, cfg.grid_west, 0, -cfg.grid_res_deg, cfg.grid_north)
    return Grid(t, w, h)


def read_multiband(path, expected=None):
    """Read a GeoTIFF into {band_name: float64 array} (nodata -> NaN).

    Band names come from the GeoTIFF descriptions when present (Earth Engine
    writes them); otherwise from `expected` (the fixed export order)."""
    path = Path(path)
    with rasterio.open(path) as src:
        data = src.read().astype("float64")
        nodata = src.nodata
        names = list(src.descriptions)
        grid = Grid(src.transform, src.width, src.height)
    if nodata is not None and not np.isnan(nodata):
        data[data == nodata] = np.nan
    data[data <= NODATA + 1] = np.nan
    if not all(names) or len(set(names)) != len(names):
        sidecar = path.with_suffix(".bands.json")
        if sidecar.exists():
            names = json.loads(sidecar.read_text())
        elif expected is not None and len(expected) == data.shape[0]:
            names = list(expected)
        else:
            raise ValueError(f"Cannot determine band names of {path}")
    return {n: data[i] for i, n in enumerate(names)}, grid


def read_stats(path):
    bands, grid = read_multiband(path, STATS_BANDS)
    missing = [b for b in STATS_BANDS if b not in bands]
    if missing:
        raise ValueError(f"{path} lacks bands {missing}")
    # Earth Engine exports masked pixels as 0 when no nodata is set; a cell
    # without observations must be NaN in every band, not 0 ppb.
    empty = ~(bands["n_obs"] > 0)
    for k in bands:
        bands[k][empty] = np.nan
    return bands, grid


def rasterize_states(ne: gpd.GeoDataFrame, grid: Grid) -> np.ndarray:
    """int array: state_id of the polygon containing each pixel centre, 0 elsewhere."""
    shapes = [(geom, int(sid)) for geom, sid in zip(ne.geometry, ne.state_id)]
    return rasterize(shapes, out_shape=(grid.height, grid.width), transform=grid.transform,
                     fill=0, all_touched=False, dtype="int16")


def land_fraction(land_geom, grid: Grid, oversample: int = 8) -> np.ndarray:
    """Share of each cell covered by a geometry (sub-pixel rasterisation)."""
    fine_t = grid.transform * Affine.scale(1.0 / oversample)
    fine = rasterize([(land_geom, 1)], out_shape=(grid.height * oversample, grid.width * oversample),
                     transform=fine_t, fill=0, all_touched=False, dtype="uint8")
    return fine.reshape(grid.height, oversample, grid.width, oversample).mean(axis=(1, 3))


def _utm_epsg(lon):
    """SIRGAS 2000 / UTM south zone EPSG code for each longitude (31982 = 22S)."""
    zone = np.floor((np.asarray(lon) + 180.0) / 6.0).astype(int) + 1
    return 31960 + zone


def distance_to_lines_km(lines: gpd.GeoDataFrame, lon, lat):
    """Distance (km) from each (lon, lat) to the nearest line and that line's index.

    Each point is measured in the SIRGAS 2000 / UTM zone that contains it
    (scale error < 0.1 % over the distances that matter here); v3 relied on a
    mixture of 1 km sampling of a 7 km image and server-side joins."""
    lon = np.asarray(lon, dtype=float).ravel()
    lat = np.asarray(lat, dtype=float).ravel()
    dist = np.full(lon.size, np.nan)
    idx = np.full(lon.size, -1, dtype=int)
    epsg = _utm_epsg(lon)
    for code in np.unique(epsg):
        sel = np.where(epsg == code)[0]
        lines_p = lines.to_crs(int(code)).geometry.values
        tree = STRtree(lines_p)
        pts = gpd.GeoSeries(points(lon[sel], lat[sel]), crs=4326).to_crs(int(code)).values
        (i_pts, i_lines), d = tree.query_nearest(pts, return_distance=True, all_matches=False)
        dist[sel[i_pts]] = d / 1000.0
        idx[sel[i_pts]] = i_lines
    return dist, idx


def distance_raster(lines: gpd.GeoDataFrame, grid: Grid, mask: np.ndarray):
    """Distance-to-nearest-line (km) and nearest-line index for pixels in `mask`."""
    d = np.full((grid.height, grid.width), np.nan)
    k = np.full((grid.height, grid.width), -1, dtype=int)
    dd, kk = distance_to_lines_km(lines, grid.lon[mask], grid.lat[mask])
    d[mask] = dd
    k[mask] = kk
    return d, k


def coast_lines_gdf(coast_geom) -> gpd.GeoDataFrame:
    g = gpd.GeoDataFrame(geometry=[coast_geom], crs=4326).explode(index_parts=False)
    return g[g.geom_type == "LineString"].reset_index(drop=True)


def write_geotiff(path, arrays: dict, grid: Grid):
    """Write named float32 bands (NaN nodata) with band descriptions."""
    names = list(arrays)
    with rasterio.open(path, "w", driver="GTiff", height=grid.height, width=grid.width,
                       count=len(names), dtype="float32", crs="EPSG:4326",
                       transform=grid.transform, nodata=np.nan, compress="deflate") as dst:
        for i, n in enumerate(names, 1):
            dst.write(np.asarray(arrays[n], dtype="float32"), i)
            dst.set_band_description(i, n)
