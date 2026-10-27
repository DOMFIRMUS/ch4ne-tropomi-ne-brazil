"""
Anomaly detection, cluster delineation and spatial attribution.

Differences from v3 (all documented in the audit report):
* the field is de-trended in time (GEE stage) and a LOCAL background is removed
  in space, so the latitudinal XCH4 gradient no longer produces "hotspots";
* the threshold uses median/MAD (robust to the retrieval-driven low tail that
  inflated sigma in v3);
* clusters are connected components with an explicit minimum size on the
  native analysis grid (v3's focal_min/focal_max on a selfMask()-ed image
  could not erode, and ran at a 10 km vectorisation scale, not the 7 km
  image grid);
* each cluster reports its PEAK pixel location and state as well as its
  centroid (v3 labelled polygon centroids as "peak" coordinates);
* distances to pipelines are real geodesic-quality distances (v3 read the
  join measure from the wrong object, so every distance was the default 0).
"""
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from rasterio.features import shapes
from scipy import ndimage, stats
from shapely.geometry import shape
from pyproj import CRS

from .raster import Grid


# ----------------------------------------------------------------- statistics
def robust_center_scale(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    return med, mad


def census_stats(x):
    """Descriptive statistics of a pixel census (population SD, as GEE's stdDev)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    med, mad = robust_center_scale(x)
    q = np.percentile(x, [2.5, 5, 25, 75, 95, 97.5])
    return dict(n=int(x.size), mean=float(x.mean()), sd=float(x.std(ddof=0)), median=float(med),
                mad_sigma=float(mad), iqr=float(q[3] - q[2]), p025=float(q[0]), p05=float(q[1]),
                p95=float(q[4]), p975=float(q[5]), skew=float(stats.skew(x)), min=float(x.min()),
                max=float(x.max()))


# ------------------------------------------------------------ local background
def window_pixels(grid: Grid, window_km: float) -> int:
    """Odd window size (pixels) spanning ~window_km in latitude."""
    px_km = grid.res * 111.195
    n = int(round(window_km / px_km))
    return n + 1 if n % 2 == 0 else n


def local_background(field, valid, grid: Grid, window_km: float):
    """Moving nan-median of `field` over `valid` pixels (square window).

    The window is expressed in pixels of constant angular size; over 1-18 deg S
    the E-W extent varies by < 5 %, which is immaterial for a smooth background."""
    arr = np.where(valid, field, np.nan)
    n = window_pixels(grid, window_km)
    return nanmedian_filter(arr, n)


def nanmedian_filter(arr, n):
    """Square moving nan-median (window n x n, NaN outside the array).

    Same result as ndimage.generic_filter(arr, np.nanmedian, size=n,
    mode='constant', cval=nan), computed row by row with a strided view:
    ~50x faster and bounded memory."""
    from numpy.lib.stride_tricks import sliding_window_view
    h = n // 2
    pad = np.pad(np.asarray(arr, float), h, mode="constant", constant_values=np.nan)
    out = np.full(arr.shape, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for r in range(arr.shape[0]):
            win = sliding_window_view(pad[r:r + n], (n, n))[0]        # (cols, n, n)
            out[r] = np.nanmedian(win.reshape(win.shape[0], -1), axis=1)
    return out


# ------------------------------------------------------------------ clusters
def label_clusters(mask, min_pixels: int, eight_connected: bool, rank_by=None):
    """Connected components >= min_pixels, relabelled 1..K by descending max(rank_by)."""
    structure = np.ones((3, 3)) if eight_connected else None
    lab, n = ndimage.label(mask, structure=structure)
    if n == 0:
        return lab, 0
    sizes = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
    keep = np.where(sizes >= min_pixels)[0] + 1
    if keep.size == 0:
        return np.zeros_like(lab), 0
    if rank_by is not None:
        peak = ndimage.maximum(np.nan_to_num(rank_by, nan=-np.inf), lab, index=keep)
        keep = keep[np.argsort(-np.asarray(peak))]
    new = np.zeros(n + 1, dtype=int)
    new[keep] = np.arange(1, keep.size + 1)
    return new[lab], int(keep.size)


def legacy_opening_clusters(mask, mode="opening", eight_connected=True):
    """Cluster counts under the possible readings of v3's morphology step
    (z.gt(2).selfMask().focal_min(1).focal_max(1), evaluated at 10 km).

    Earth Engine reducers ignore masked pixels, so on a selfMask()-ed image
    focal_min can never erode: the step was either the identity or a
    dilation by two pixels (5 x 5), never the "closing" (or opening) the
    manuscript describes. mode: 'identity', 'dilate2' or 'opening' (intended)."""
    st = np.ones((3, 3), bool)
    if mode == "opening":
        work = ndimage.binary_opening(mask, structure=st)
    elif mode == "dilate2":
        work = ndimage.binary_dilation(mask, structure=st, iterations=2)
    else:
        work = mask
    lab, n = ndimage.label(work, structure=np.ones((3, 3)) if eight_connected else None)
    return lab, n


def polygonize(labels, grid: Grid) -> gpd.GeoDataFrame:
    geoms, ids = [], []
    for geom, val in shapes(labels.astype("int32"), mask=labels > 0, transform=grid.transform):
        geoms.append(shape(geom))
        ids.append(int(val))
    gdf = gpd.GeoDataFrame({"cluster": ids}, geometry=geoms, crs=4326)
    return gdf.dissolve(by="cluster").reset_index()


def _aeqd(lon, lat):
    return CRS.from_proj4(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m")


def cluster_table(labels, fields: dict, grid: Grid, state_raster, ne: gpd.GeoDataFrame,
                  lines: gpd.GeoDataFrame, coast_dist_km, pipe_dist_km, polys: gpd.GeoDataFrame):
    """One row per cluster with peak/centroid positions, state, sizes and distances.

    fields must contain 'z' (detection statistic), 'enh' (ppb enhancement over
    local background), 'xch4_mean' and 'n_obs'."""
    area = grid.pixel_area_km2
    lon, lat = grid.lon, grid.lat
    uf_by_id = dict(zip(ne.state_id, ne.UF))
    rows = []
    for k in range(1, int(labels.max()) + 1):
        m = labels == k
        if not m.any():
            continue
        z = np.where(m, fields["z"], -np.inf)
        r, c = np.unravel_index(np.argmax(z), z.shape)
        w = area[m]
        poly = polys.loc[polys.cluster == k, "geometry"].iloc[0]
        cen_lon = float(np.sum(lon[m] * w) / w.sum())
        cen_lat = float(np.sum(lat[m] * w) / w.sum())
        # geodesic-quality polygon-to-line distance in an azimuthal equidistant
        # projection centred on the cluster
        aeqd = _aeqd(cen_lon, cen_lat)
        poly_p = gpd.GeoSeries([poly], crs=4326).to_crs(aeqd).iloc[0]
        lines_p = lines.to_crs(aeqd)
        d_all = lines_p.distance(poly_p)
        j = int(d_all.idxmin())
        ids_in = np.unique(state_raster[m])
        ids_in = ids_in[ids_in > 0]
        rows.append(dict(
            cluster=k,
            n_pixels=int(m.sum()),
            area_km2=float(w.sum()),
            peak_lat=float(lat[r, c]), peak_lon=float(lon[r, c]),
            peak_UF=uf_by_id.get(int(state_raster[r, c]), "outside NE"),
            UFs_spanned=",".join(sorted(uf_by_id[int(i)] for i in ids_in)),
            centroid_lat=cen_lat, centroid_lon=cen_lon,
            peak_z=float(fields["z"][r, c]),
            peak_enh_ppb=float(fields["enh"][r, c]),
            mean_enh_ppb=float(np.nanmean(fields["enh"][m])),
            peak_xch4_ppb=float(fields["xch4_mean"][r, c]),
            peak_n_obs=float(fields["n_obs"][r, c]),
            min_n_obs=float(np.nanmin(fields["n_obs"][m])),
            dist_poly_to_pipeline_km=float(d_all.min() / 1000.0),
            dist_peak_to_pipeline_km=float(pipe_dist_km[r, c]),
            nearest_asset=str(lines.loc[j, "asset_name"]),
            nearest_asset_class=str(lines.loc[j, "asset_class"]),
            dist_peak_to_coast_km=float(coast_dist_km[r, c]),
        ))
    return pd.DataFrame(rows)


def mirror_qvalues(peak_z_pos, peak_z_neg):
    """Per-cluster false-discovery q-values from the mirror (negative) clusters.

    For a positive cluster with peak z = t:
        q(t) = min over t' <= t of (1 + #negative peaks with |z| >= t') / #positive peaks >= t'
    i.e. the knockoff+ / target-decoy "+1" estimate (Barber & Candes 2015,
    Ann. Stat. 43, 2055); without the +1 the top clusters would get q = 0
    whatever the evidence. Consequence: q <= 0.10 needs at least 10 positive
    clusters at or above t, so with fewer than 10 detections none can be
    declared significant by construction.
    This is a HEURISTIC false-discovery estimate, not a guaranteed FDR bound:
    it is calibrated only if, in the absence of emission, exceedances above
    and below the local background are equally likely everywhere (sign
    symmetry). Retrieval artefacts need not be symmetric and their sign can
    change with data density or distance to the coast (checked per stratum in
    TableS14/TableS15), so q can be conservative in some strata and liberal in
    others."""
    pos = np.sort(np.asarray(peak_z_pos, float))[::-1]
    neg = np.sort(np.abs(np.asarray(peak_z_neg, float)))
    q = np.array([min(1.0, (1 + neg.size - np.searchsorted(neg, t, side="left")) / (i + 1))
                  for i, t in enumerate(pos)])
    q = np.minimum.accumulate(q[::-1])[::-1]
    order = np.argsort(-np.asarray(peak_z_pos, float), kind="stable")
    out = np.empty_like(q)
    out[order] = q
    return out


def cluster_peaks(labels, z):
    k = int(labels.max())
    if k == 0:
        return np.array([])
    return np.asarray(ndimage.maximum(np.nan_to_num(z, nan=-np.inf), labels, index=np.arange(1, k + 1)))


def cluster_peak_values_at(labels, score, field):
    """For each label 1..K, the value of `field` at the pixel where `score` peaks."""
    k = int(labels.max())
    if k == 0:
        return np.array([])
    pos = ndimage.maximum_position(np.nan_to_num(score, nan=-np.inf), labels, index=np.arange(1, k + 1))
    return np.array([field[r, c] for r, c in pos], dtype=float)


# --------------------------------------------------------------- attribution
def attribution_enrichment(peak_dist_pos, peak_dist_neg, eligible, pipe_dist_km, area_km2, radii,
                           peak_coast_pos=None, peak_coast_neg=None, coast_strata_km=(50.0,)):
    """How often cluster peaks lie within r km of a pipeline, against two nulls.

    1. Area share (descriptive base rate): under spatial randomness a peak
       falls within r km with probability p_r = eligible area within r km /
       eligible area; the count is Binomial(N, p_r). This null ignores that
       detectability varies in space (noisier, coastal cells cross the
       threshold more often), so it is optimistic.
    2. Mirror clusters (empirical null): the negative clusters obtained with
       identical rules. A one-sided Fisher exact test compares the share of
       positive and of mirror peaks within r km. Because mirror and positive
       peaks can be distributed differently with respect to the coast (where
       pipelines run), the comparison is also STRATIFIED by the peak's
       distance to the coast (strata split at coast_strata_km): Mantel-Haenszel
       pooled odds ratio and Cochran-Mantel-Haenszel test (two-sided), with the
       per-stratum counts. Strata without positive or without mirror peaks
       carry no information and are left out of the pooled estimate.
    Peak distances are pixel-centre distances, matching both nulls."""
    from statsmodels.stats.contingency_tables import StratifiedTable
    pos = np.asarray(peak_dist_pos, float)
    neg = np.asarray(peak_dist_neg, float)
    N, M = pos.size, neg.size
    tot = np.nansum(np.where(eligible, area_km2, 0))
    edges = [-np.inf] + list(coast_strata_km) + [np.inf]
    cp = None if peak_coast_pos is None else np.asarray(peak_coast_pos, float)
    cn = None if peak_coast_neg is None else np.asarray(peak_coast_neg, float)
    out = []
    for r in radii:
        near_area = np.nansum(np.where(eligible & (pipe_dist_km <= r), area_km2, 0))
        p = near_area / tot if tot > 0 else np.nan
        k = int(np.sum(pos <= r))
        k_neg = int(np.sum(neg <= r))
        exp = N * p
        p_binom = float(stats.binom.sf(k - 1, N, p)) if N and 0 < p < 1 else np.nan
        if N and M:
            odds, p_fisher = stats.fisher_exact([[k, N - k], [k_neg, M - k_neg]], alternative="greater")
        else:
            odds, p_fisher = np.nan, np.nan
        row = dict(radius_km=r, n_clusters=N, observed_peaks_within=k, share_within=k / N if N else np.nan,
                   area_share_within=p, expected_peaks_within_area_null=exp,
                   enrichment_vs_area=(k / exp) if exp > 0 else np.nan, p_binomial_area_null=p_binom,
                   n_mirror_clusters=M, mirror_peaks_within=k_neg,
                   mirror_share_within=k_neg / M if M else np.nan,
                   odds_ratio_vs_mirror=float(odds), p_fisher_vs_mirror=float(p_fisher))
        if cp is not None and cn is not None:
            tabs, desc = [], []
            for lo, hi in zip(edges[:-1], edges[1:]):
                sp, sn = (cp > lo) & (cp <= hi), (cn > lo) & (cn <= hi)
                a_, b_ = int(np.sum(pos[sp] <= r)), int(sp.sum())
                c_, d_ = int(np.sum(neg[sn] <= r)), int(sn.sum())
                lab = f"coast {'' if not np.isfinite(lo) else f'>{lo:g}'}{'' if not np.isfinite(hi) else f'<={hi:g}'} km"
                desc.append(f"{lab}: pos {a_}/{b_}, mirror {c_}/{d_}")
                t_ = np.array([[a_, b_ - a_], [c_, d_ - c_]], float)
                # informative only with both groups present and peaks both within and beyond r
                if b_ and d_ and t_[:, 0].sum() > 0 and t_[:, 1].sum() > 0:
                    tabs.append(t_)
            if len(tabs) == 1:   # one informative stratum: exact test (Haldane OR if a cell is 0)
                t_ = tabs[0]
                tt = t_ + 0.5 if (t_ == 0).any() else t_
                or_mh = float(tt[0, 0] * tt[1, 1] / (tt[0, 1] * tt[1, 0]))
                p_cmh = float(stats.fisher_exact(t_.astype(int), alternative="two-sided")[1])
                test = "fisher_exact_single_informative_stratum"
            elif len(tabs) > 1:
                st = StratifiedTable(tabs)
                or_mh = float(st.oddsratio_pooled)
                p_cmh = float(st.test_null_odds(correction=True).pvalue)
                test = "cochran_mantel_haenszel"
            else:
                or_mh, p_cmh, test = np.nan, np.nan, "no_informative_stratum"
            row.update(odds_ratio_mh_coast_stratified=or_mh, p_cmh_coast_stratified=p_cmh,
                       coast_stratified_test=test, coast_strata_counts="; ".join(desc))
        out.append(row)
    return pd.DataFrame(out)
