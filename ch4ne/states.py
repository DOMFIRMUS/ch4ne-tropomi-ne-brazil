"""
Per-state statistics (manuscript Table 2 / violin figure) with the state
identity taken from the polygon attributes, plus the independent cross-check
suggested in PENDENCIAS: every cluster peak must fit under the maximum of the
state that contains it, where "contains" is decided by a vector
point-in-polygon test, independent of the raster path that built the table.
"""
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point


def state_table(fields: dict, state_raster, ne: gpd.GeoDataFrame, valid, grid, regional_mean_of: dict):
    """One row per state. `regional_mean_of[var]` is the NE-wide census mean used
    for the 'anomaly vs region' columns (same scale, same mask)."""
    rows = []
    lon, lat = grid.lon, grid.lat
    for _, s in ne.iterrows():
        m = (state_raster == s.state_id) & valid
        if not (m & np.isfinite(fields["xch4_mean"])).any():   # e.g. a small state in a short robustness period
            rows.append(dict(UF=s.UF, state_name=s.state_name, n_pixels=0))
            continue
        x = fields["xch4_mean"][m]
        a = fields["anom_mean"][m]
        e = fields["enh"][m]
        i_max = np.nanargmax(np.where(m, fields["xch4_mean"], -np.inf))
        r, c = np.unravel_index(i_max, m.shape)
        rows.append(dict(
            UF=s.UF, state_name=s.state_name, n_pixels=int(m.sum()),
            xch4_mean_ppb=np.nanmean(x), xch4_median_ppb=np.nanmedian(x), xch4_sd_ppb=np.nanstd(x),
            xch4_max_ppb=np.nanmax(x), xch4_max_lat=float(lat[r, c]), xch4_max_lon=float(lon[r, c]),
            xch4_minus_region_ppb=np.nanmean(x) - regional_mean_of["xch4_mean"],
            anom_mean_ppb=np.nanmean(a), anom_minus_region_ppb=np.nanmean(a) - regional_mean_of["anom_mean"],
            enh_median_ppb=np.nanmedian(e), enh_sd_ppb=np.nanstd(e),
            temporal_sd_median_ppb=np.nanmedian(fields["anom_sd"][m]),
            n_obs_median=np.nanmedian(fields["n_obs"][m]),
        ))
    return pd.DataFrame(rows).sort_values("xch4_mean_ppb", ascending=False).reset_index(drop=True)


def check_peaks_within_state_max(clusters: pd.DataFrame, states_tab: pd.DataFrame,
                                 ne: gpd.GeoDataFrame, peak_col="peak_xch4_ppb", tol=1e-6):
    """Vector point-in-polygon cross-check. Returns a list of violations (empty = pass)."""
    mx = states_tab.set_index("UF")["xch4_max_ppb"]
    bad = []
    for _, r in clusters.iterrows():
        hit = ne[ne.contains(Point(r.peak_lon, r.peak_lat))]
        if hit.empty or not np.isfinite(r[peak_col]):
            continue
        uf = hit.iloc[0].UF
        if r[peak_col] > mx[uf] + tol:
            bad.append(dict(cluster=int(r.cluster), UF=uf, peak=float(r[peak_col]), state_max=float(mx[uf])))
    return bad
