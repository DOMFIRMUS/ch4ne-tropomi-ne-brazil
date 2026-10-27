"""
Local analysis stage: reads the rasters exported by gee_export.py and writes
every table, figure and number used in the manuscript, plus automatic
consistency checks. It refuses to run on anything that is not a verified
Earth Engine export (see provenance.py).
"""
import importlib.metadata as md
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from . import figures as F
from .boundaries import coastline, load_states, ne_states
from .config import Config
from .detection import (attribution_enrichment, census_stats, cluster_peak_values_at, cluster_peaks, cluster_table,
                        label_clusters, legacy_opening_clusters, local_background, mirror_qvalues, polygonize,
                        robust_center_scale)
from .infra import inventory_table, load_infrastructure, onshore_lines
from .provenance import verify as verify_provenance
from .raster import (coast_lines_gdf, distance_raster, land_fraction, rasterize_states, read_multiband,
                     read_stats, write_geotiff)
from .states import check_peaks_within_state_max, state_table
from .timeseries import cluster_monthly_series, fit_trend, read_background_csv
from .zones import (COVARIATES, assign_rings, block_bootstrap_compare, block_ids, coast_collocation,
                    leave_one_block_out, matched_contrast, matched_controls, regressions, ring_table,
                    stratified_contrast)

PROX_VARS = ["xch4_mean", "bg_sampled", "anom_mean", "enh", "anom_sd", "tail_up_90", "tail_lo_10", "tail_asym_90",
             "tail_up_95", "tail_lo_05", "n_obs", "mean_t_years", "dist_coast_km", "albedo_b7", "water_frac"]
# Inferential contrasts use de-trended variables only. Raw xch4_mean mixes in
# the sampling time of each cell (growth ~8-10 ppb/yr); its sampling part,
# bg_sampled = xch4_mean - anom_mean, is reported as a diagnostic row.
CONTRAST_VARS = ["anom_mean", "enh", "anom_sd", "tail_up_90", "tail_lo_10", "tail_asym_90"]
DIAGNOSTIC_VARS = ["bg_sampled"]
REG_VARS = ["anom_mean", "enh", "anom_sd", "tail_up_90", "tail_lo_10", "tail_asym_90"]
SOFTWARE = ("numpy", "pandas", "scipy", "geopandas", "shapely", "pyproj", "pyogrio", "rasterio", "affine",
            "statsmodels", "matplotlib", "earthengine-api")


def _json_safe(o):
    """Strict JSON: non-finite floats become null, numpy scalars become Python values."""
    if isinstance(o, dict):
        return {str(k): _json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_safe(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (str, int, bool)) or o is None:
        return o
    return str(o)


def _csv(df, path):
    df.to_csv(path, index=False, encoding="utf-8", float_format="%.6g")
    return path


def _ver(pkg):
    try:
        return md.version(pkg)
    except md.PackageNotFoundError:
        return None


def run_record(cfg: Config, extra_cli=None):
    import pyproj
    import rasterio
    software = {p: _ver(p) for p in SOFTWARE}
    software.update(python=sys.version, platform=platform.platform(), ch4ne=__version__,
                    gdal=rasterio.__gdal_version__, proj=pyproj.proj_version_str)
    return dict(created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                config=cfg.as_dict(), cli=extra_cli if extra_cli is not None else sys.argv, software=software)


def clean_outputs(out: Path):
    """Remove this run's previous products (never out/gee), so conditional
    outputs of an earlier run cannot mix with the new ones."""
    for pat in ("tables/*.csv", "figures/Fig*.png", "figures/Fig*.pdf"):
        for p in out.glob(pat):
            p.unlink()
    for name in ("clusters_v4.geojson", "manuscript_numbers.json", "checks.json", "derived_fields.tif",
                 "pixels_strict.csv.gz"):
        p = out / name
        if p.exists():
            p.unlink()


def months_with_data(gee_dir, grid):
    """Number of months with at least one retrieval per cell, recounted from
    the monthly count stack (cross-check of the n_months band of stats.tif)."""
    MC, g = read_multiband(Path(gee_dir) / "monthly_count.tif")
    if (g.width, g.height) != (grid.width, grid.height):
        raise RuntimeError("monthly_count.tif is not on the analysis grid")
    return np.sum(np.stack([np.nan_to_num(v) > 0 for v in MC.values()]), axis=0).astype(float)


def detect(f, grid, screen, bg_screen, cfg, z_threshold=None, min_obs=None, min_months=None, window_km=None):
    """Local-background removal, robust z, connected components, mirror q-values.

    `screen` / `bg_screen` are surface/region masks WITHOUT the data-count
    rules, which are applied here so the sensitivity runs can vary them."""
    z_thr = cfg.z_threshold if z_threshold is None else z_threshold
    mo = cfg.min_obs if min_obs is None else min_obs
    mm = cfg.min_valid_months if min_months is None else min_months
    enough = f["n_obs"] >= mo
    if mm > 0:
        enough &= f["n_months"] >= mm
    bgv = bg_screen & enough
    st = screen & enough
    bg = local_background(f["anom_mean"], bgv, grid, window_km or cfg.local_bg_window_km)
    enh = f["anom_mean"] - bg
    med, mad = robust_center_scale(enh[st])
    z = (enh - med) / mad
    det = st & (z > z_thr)
    labels, k = label_clusters(det, cfg.min_cluster_pixels, cfg.eight_connected, rank_by=np.where(st, z, np.nan))
    # Mirror test: identical rules on z < -threshold. Under noise symmetric
    # about the local background, mirror clusters estimate false discoveries.
    neg_labels, k_neg = label_clusters(st & (z < -z_thr), cfg.min_cluster_pixels, cfg.eight_connected)
    q = mirror_qvalues(cluster_peaks(labels, z), cluster_peaks(neg_labels, -z)) if k else np.array([])
    return dict(bg=bg, enh=enh, z=z, med=med, mad=mad, det=det, labels=labels, n_clusters=k,
                neg_labels=neg_labels, n_negative_mirror=k_neg, q=q,
                n_q=int(np.sum(q <= cfg.q_threshold)), screen=st, bg_screen=bgv)


def peak_distances(D, pipe_d):
    """Pipeline distance at the peak of every positive and every mirror cluster."""
    pos = cluster_peak_values_at(D["labels"], D["z"], pipe_d)
    neg = cluster_peak_values_at(D["neg_labels"], -D["z"], pipe_d)
    return pos, neg


def mirror_strata(D, coast_d, n_obs, cfg, coast_km=50.0, nobs_split=None):
    """Mirror q-values recomputed within strata of the peak's distance to the
    coast and data density (split at the median n_obs of the screened cells
    unless given): a check of the sign-symmetry assumption."""
    if nobs_split is None:
        nobs_split = float(np.nanmedian(n_obs[D["screen"]]))
    zp, zn = cluster_peaks(D["labels"], D["z"]), cluster_peaks(D["neg_labels"], -D["z"])
    cp, cn = peak_distances(D, coast_d)
    npk = cluster_peak_values_at(D["labels"], D["z"], n_obs)
    nnk = cluster_peak_values_at(D["neg_labels"], -D["z"], n_obs)
    fixed = 3.0 * cfg.min_obs        # a fixed split as well, to show how much q depends on the cut
    strata = [("all", np.ones(zp.size, bool), np.ones(zn.size, bool)),
              (f"peak <= {coast_km:g} km from coast", cp <= coast_km, cn <= coast_km),
              (f"peak > {coast_km:g} km from coast", cp > coast_km, cn > coast_km),
              (f"peak n_obs < {nobs_split:g} (median)", npk < nobs_split, nnk < nobs_split),
              (f"peak n_obs >= {nobs_split:g} (median)", npk >= nobs_split, nnk >= nobs_split),
              (f"peak n_obs >= {fixed:g} (3 x min_obs)", npk >= fixed, nnk >= fixed)]
    rows = []
    for name, sp, sn in strata:
        q = mirror_qvalues(zp[sp], zn[sn]) if sp.any() else np.array([])
        rows.append(dict(stratum=name, n_positive=int(sp.sum()), n_mirror=int(sn.sum()),
                         fdr_point_estimate=(sn.sum() / sp.sum()) if sp.any() else np.nan,
                         min_q=float(q.min()) if q.size else np.nan,
                         n_q_le_threshold=int(np.sum(q <= cfg.q_threshold))))
    return pd.DataFrame(rows)


def sign_symmetry(D, coast_d, n_obs, cfg):
    """Share of screened cells above +z and below -z, by data density and by
    distance to the coast. Under the mirror's assumption both shares match."""
    st, z, thr = D["screen"], D["z"], cfg.z_threshold
    q = np.nanquantile(n_obs[st], [0.2, 0.4, 0.6, 0.8])      # quintiles of THIS run's screened cells
    nobs_edges = tuple([float(np.nanmin(n_obs[st]))] + [float(x) for x in q] + [np.inf])
    rows = []
    for var, edges, vals in (("n_obs", nobs_edges, n_obs),
                             ("dist_coast_km", (0, 20, 50, 150, np.inf), coast_d)):
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = st & (vals >= lo) & (vals < hi)
            n = int(m.sum())
            if n:
                zm = z[m]
                med = float(np.median(zm))
                up, dn = float(np.mean(zm > thr)), float(np.mean(zm < -thr))
                upc, dnc = float(np.mean(zm - med > thr)), float(np.mean(zm - med < -thr))
            else:
                med = up = dn = upc = dnc = np.nan
            rows.append(dict(variable=var, low=lo, high=hi, n_cells=n, median_z=med,
                             share_z_above=up, share_z_below=dn,
                             ratio_above_to_below=up / dn if dn else np.nan,
                             ratio_after_centring_on_stratum_median=upc / dnc if dnc else np.nan))
    return pd.DataFrame(rows)


def significant_clusters(d, grid, pipe_d, coast_d, n_obs, cfg):
    """Peak position, q and distances of the clusters with q <= threshold."""
    zp = cluster_peaks(d["labels"], d["z"])
    keep = np.where(d["q"] <= cfg.q_threshold)[0]
    lon = cluster_peak_values_at(d["labels"], d["z"], grid.lon)
    lat = cluster_peak_values_at(d["labels"], d["z"], grid.lat)
    pdist = cluster_peak_values_at(d["labels"], d["z"], pipe_d)
    cdist = cluster_peak_values_at(d["labels"], d["z"], coast_d)
    nob = cluster_peak_values_at(d["labels"], d["z"], n_obs)
    return pd.DataFrame(dict(cluster=keep + 1, peak_lon=lon[keep], peak_lat=lat[keep], peak_z=zp[keep],
                             q_mirror=d["q"][keep], dist_peak_to_pipeline_km=pdist[keep],
                             dist_peak_to_coast_km=cdist[keep], peak_n_obs=nob[keep]))


def pixel_frame(f, grid, mask, extra: dict):
    d = dict(lon=grid.lon[mask], lat=grid.lat[mask])
    for k, v in f.items():
        d[k] = v[mask]
    for k, v in extra.items():
        d[k] = v[mask]
    df = pd.DataFrame(d)
    df["tail_asym_90"] = df["tail_up_90"] - df["tail_lo_10"]   # > 0: upper tail wider (intermittent emission)
    df["log_n_obs"] = np.log(df["n_obs"])
    df["bg_sampled"] = df["xch4_mean"] - df["anom_mean"]   # background at this cell's sampling times
    return df


def legacy_reproduction(lx, lgrid, ne, cfg):
    """v3 statistics recomputed with Earth Engine's reduceRegion inclusion rule
    (weighted reducers: cells with >= 1/256 of their area in the region,
    weighted by the covered fraction), next to the plain pixel-centre census."""
    ne_union = ne.union_all()
    fr = land_fraction(ne_union, lgrid, oversample=16)
    centroid = (rasterize_states(ne, lgrid) > 0) & np.isfinite(lx)
    touched = (fr >= 1 / 256) & np.isfinite(lx)
    w = np.where(touched, fr, 0.0)
    mu_w = float(np.sum(w * np.nan_to_num(lx)) / np.sum(w))
    sd_w = float(np.sqrt(np.sum(w * (np.nan_to_num(lx) - mu_w) ** 2) / np.sum(w)))
    cs = census_stats(lx[centroid])
    return dict(mu_weighted=mu_w, sigma_weighted=sd_w, n_centroid=int(centroid.sum()), n_touched=int(touched.sum()),
                n_weighted_sum=float(w.sum()), mu_centroid=cs["mean"], sigma_centroid=cs["sd"]), centroid


def run(cfg: Config, gee_dir, out_dir=None, sensitivity=True, cli=None):
    gee_dir = Path(gee_dir)
    provenance = verify_provenance(gee_dir, cfg)  # refuses anything not downloaded from Earth Engine
    out = Path(out_dir) if out_dir else cfg.out_dir
    T, FG = out / "tables", out / "figures"
    T.mkdir(parents=True, exist_ok=True)
    FG.mkdir(parents=True, exist_ok=True)
    clean_outputs(out)
    F.setup_style(cfg)
    numbers = {"run": run_record(cfg, cli)}
    checks = {}

    # ------------------------------------------------ boundaries & infrastructure
    states = load_states(cfg)
    ne = ne_states(states, cfg)
    infra = load_infrastructure(cfg, states, ne)
    lines = onshore_lines(infra, states, cfg)
    _csv(inventory_table(infra), T / "Table1_infrastructure_inventory.csv")
    infra.drop(columns="geometry").to_csv(T / "TableS1_infrastructure_segments.csv", index=False)
    numbers["infra_segments_total_selected"] = int(len(infra))
    numbers["infra_segments_by_location"] = infra["location"].value_counts().to_dict()
    numbers["infra_length_on_ne_land_km_geodesic"] = float(infra["length_ne_land_km"].sum())

    # ----------------------------------------------------------------- rasters
    f, grid = read_stats(gee_dir / "stats.tif")
    state_r = rasterize_states(ne, grid)
    in_ne = state_r > 0
    land = land_fraction(states.union_all(), grid)
    has = np.isfinite(f["xch4_mean"])
    landish = in_ne | (land >= 0.5)
    pipe_d, _ = distance_raster(lines, grid, in_ne)
    coast_d, _ = distance_raster(coast_lines_gdf(coastline(states)), grid, landish)  # NE and neighbours

    # months with data: computed on Earth Engine in the statistics pass;
    # recounted from the monthly stack when it was exported
    has_monthly = (gee_dir / "monthly_count.tif").exists() and (gee_dir / "monthly_anom.tif").exists()
    if has_monthly:
        nm = months_with_data(gee_dir, grid)
        checks["n_months_equals_monthly_stack"] = bool(np.array_equal(f["n_months"][has], nm[has]))
    obs_ok = has & (f["n_obs"] >= cfg.min_obs) & (f["n_months"] >= cfg.min_valid_months)
    surface_ok = (land >= cfg.min_land_fraction) & ~(f["water_frac"] > cfg.max_water_fraction) \
        & (np.nan_to_num(coast_d, nan=0.0) >= cfg.coast_buffer_km)
    # region/surface screens without the data-count rules; the SAME surface
    # screen is used for the local background, so coastal cells cannot pull
    # the background of near-pipeline cells down
    relaxed_s = has & in_ne
    strict_s = relaxed_s & surface_ok
    bg_s = has & surface_ok
    relaxed = relaxed_s & obs_ok
    strict = strict_s & obs_ok
    numbers.update(
        min_valid_months_applied=cfg.min_valid_months,
        pixels_in_NE_with_data=int((has & in_ne).sum()),
        pixels_relaxed_screen=int(relaxed.sum()), pixels_strict_screen=int(strict.sum()),
        pixels_removed_low_counts=int((has & in_ne & ~obs_ok).sum()),
        pixels_removed_coast_water=int((relaxed & ~strict).sum()),
        n_obs_median_NE=float(np.nanmedian(f["n_obs"][has & in_ne])),
        n_obs_p05_NE=float(np.nanpercentile(f["n_obs"][has & in_ne], 5)),
        n_months_median_NE=float(np.nanmedian(f["n_months"][has & in_ne])),
    )

    # ------------------------------------------------ legacy (v3) reproduction
    leg_path = gee_dir / "legacy.tif"
    legacy_x, mu, sd = f["xch4_mean"][relaxed], np.nan, np.nan
    if leg_path.exists() and "legacy" in provenance.get("files", {}):   # never a stale, unverified file
        L, lgrid = read_multiband(leg_path, ["legacy_xch4"])
        lx = L["legacy_xch4"]
        rep, lin = legacy_reproduction(lx, lgrid, ne, cfg)
        mu, sd = rep["mu_weighted"], rep["sigma_weighted"]
        lmask = lin & ((lx - mu) / sd > 2.0)
        n_by_mode = {m: legacy_opening_clusters(lmask, mode=m)[1] for m in ("identity", "dilate2", "opening")}
        leg = dict(rep, threshold=mu + 2 * sd, area_km2=float(lgrid.pixel_area_km2[lmask].sum()),
                   n_pixels_z2=int(lmask.sum()), clusters_identity=n_by_mode["identity"],
                   clusters_dilate2=n_by_mode["dilate2"], clusters_opening=n_by_mode["opening"],
                   published_mu=cfg.published_mu, published_sigma=cfg.published_sigma, published_n=cfg.published_n)
        pd.DataFrame([leg]).to_csv(T / "TableS2_legacy_v3_reproduction.csv", index=False)
        numbers["legacy"] = leg
        checks["legacy_reproduces_published_mu_sigma"] = bool(
            abs(mu - cfg.published_mu) < 0.05 and abs(sd - cfg.published_sigma) < 0.05)
        checks["legacy_published_n_between_centroid_and_touched"] = bool(
            rep["n_centroid"] <= cfg.published_n <= rep["n_touched"])
        legacy_x = lx[lin]

    # ------------------------------------------------------------- detection
    D = detect(f, grid, strict_s, bg_s, cfg)
    labels, enh, z = D["labels"], D["enh"], D["z"]
    polys = polygonize(labels, grid) if D["n_clusters"] else None
    fields = dict(z=z, enh=enh, xch4_mean=f["xch4_mean"], n_obs=f["n_obs"])
    ct = cluster_table(labels, fields, grid, state_r, ne, lines, coast_d, pipe_d, polys) if polys is not None \
        else pd.DataFrame()
    if len(ct):
        ct["q_mirror"] = D["q"][ct["cluster"].values - 1]
        ct[f"significant_q{int(round(cfg.q_threshold * 100)):02d}"] = ct["q_mirror"] <= cfg.q_threshold
    _csv(ct, T / "Table3_clusters.csv")
    if polys is not None:
        polys.merge(ct, on="cluster").to_file(out / "clusters_v4.geojson", driver="GeoJSON")
    pos_d, neg_d = peak_distances(D, pipe_d)
    pos_c, neg_c = peak_distances(D, coast_d)
    enrich = attribution_enrichment(pos_d, neg_d, strict, pipe_d, grid.pixel_area_km2, cfg.attribution_radius_km,
                                    peak_coast_pos=pos_c, peak_coast_neg=neg_c,
                                    coast_strata_km=cfg.attribution_coast_strata_km)
    _csv(enrich, T / "Table4_attribution_enrichment.csv")
    _csv(mirror_strata(D, coast_d, f["n_obs"], cfg), T / "TableS14_mirror_by_stratum.csv")
    _csv(sign_symmetry(D, coast_d, f["n_obs"], cfg), T / "TableS15_sign_symmetry.csv")
    area = grid.pixel_area_km2
    numbers.update(
        enh_median_ppb=D["med"], enh_robust_sigma_ppb=D["mad"],
        detection_threshold_enh_ppb=D["med"] + cfg.z_threshold * D["mad"],
        n_clusters=D["n_clusters"], n_negative_mirror_clusters=D["n_negative_mirror"],
        fdr_mirror_point_estimate=(D["n_negative_mirror"] / D["n_clusters"]) if D["n_clusters"] else None,
        n_clusters_q_le_threshold=D["n_q"], q_threshold=cfg.q_threshold,
        area_suprathreshold_km2=float(area[D["det"]].sum()),
        area_in_clusters_km2=float(area[labels > 0].sum()),
        area_strict_km2=float(area[strict].sum()),
    )
    write_geotiff(out / "derived_fields.tif",
                  dict(enh=enh, z=z, local_bg=D["bg"], cluster=labels.astype(float),
                       mirror_cluster=D["neg_labels"].astype(float), dist_pipeline_km=pipe_d,
                       dist_coast_km=coast_d, land_frac=land, n_months=f["n_months"],
                       screen_strict=strict.astype(float), state_id=state_r.astype(float)), grid)

    # -------------------------------------------------------------- states
    extra = dict(enh=enh, z=z, dist_pipe_km=pipe_d, dist_coast_km=coast_d, land_frac=land,
                 state_id=state_r.astype(float))
    uf = dict(zip(ne.state_id, ne.UF))
    tabs = {}
    for name, mask in (("strict", strict), ("relaxed", relaxed)):
        reg = {"xch4_mean": float(np.nanmean(f["xch4_mean"][mask])),
               "anom_mean": float(np.nanmean(f["anom_mean"][mask]))}
        tabs[name] = state_table(dict(f, enh=enh), state_r, ne, mask, grid, reg)
    _csv(tabs["strict"], T / "Table2_states_strict.csv")
    _csv(tabs["relaxed"], T / "TableS3_states_relaxed.csv")
    viol = check_peaks_within_state_max(ct, tabs["strict"], ne) if len(ct) else []
    checks["cluster_peaks_within_state_max"] = not viol
    checks["cluster_peak_violations"] = viol

    # ----------------------------------------------------------- proximity
    colloc = coast_collocation(pipe_d, coast_d, in_ne)
    _csv(colloc, T / "TableS9_pipeline_coast_collocation.csv")
    numbers["pipeline_coast_collocation"] = colloc.to_dict(orient="records")
    results = {}
    for name, mask in (("strict", strict), ("relaxed", relaxed)):
        df = pixel_frame(f, grid, mask, extra)
        df["UF"] = df["state_id"].astype(int).map(uf)
        rings, order = assign_rings(df["dist_pipe_km"].values, cfg.ring_edges_km)
        df["ring"] = rings
        df["block"] = block_ids(df["lon"], df["lat"], cfg.bootstrap_block_deg)
        checks[f"rings_exclusive_and_complete_{name}"] = bool(pd.notna(df["ring"]).all())
        rt = ring_table(df, PROX_VARS, order=order)
        _csv(rt, T / f"Table5_proximity_rings_{name}.csv")
        comps = pd.concat([block_bootstrap_compare(df, v, cfg.near_ring, cfg.far_ring, cfg.n_boot, cfg.random_seed)
                           .assign(role="diagnostic" if v in DIAGNOSTIC_VARS else "test")
                           for v in CONTRAST_VARS + DIAGNOSTIC_VARS], ignore_index=True)
        _csv(comps, T / f"Table6_near_vs_far_bootstrap_{name}.csv")
        results[name] = dict(df=df, order=order, rings=rt)

    df = results["strict"]["df"]
    covs = [c for c in COVARIATES if df[c].notna().mean() > 0.95]
    dm = df.dropna(subset=covs).reset_index(drop=True)
    near_m = (dm["ring"] == cfg.near_ring).values
    far_m = (dm["ring"] == cfg.far_ring).values
    extra = ("mean_t_years", "bg_sampled")
    cal = cfg.match_caliper if cfg.match_method == "euclidean" else cfg.match_caliper_ps_sd
    ctrl, bal, matched, msum = matched_controls(dm, near_m, far_m, covs, cfg.n_matched_controls,
                                                cal, balance_extra=extra, method=cfg.match_method)
    _csv(bal, T / "TableS4_matching_balance.csv")
    numbers["matching"] = msum
    mc = pd.concat([matched_contrast(dm, ctrl, matched, v, cfg.n_boot, cfg.random_seed, covariates=covs,
                                     far_mask=far_m) for v in CONTRAST_VARS], ignore_index=True)
    _csv(mc, T / "Table7_near_vs_matched_controls.csv")
    # sensitivity: Euclidean matching in the standardised covariate space
    ctrl_e, bal_e, matched_e, msum_e = matched_controls(dm, near_m, far_m, covs, cfg.n_matched_controls,
                                                        cfg.match_caliper, balance_extra=extra, method="euclidean")
    numbers["matching_euclidean_sensitivity"] = msum_e
    mce = pd.concat([matched_contrast(dm, ctrl_e, matched_e, v, cfg.n_boot, cfg.random_seed, covariates=covs,
                                      far_mask=far_m) for v in CONTRAST_VARS], ignore_index=True)
    _csv(pd.concat([mce.assign(method="euclidean"), bal_e.assign(method="euclidean")], ignore_index=True),
         T / "TableS13_matching_euclidean_sensitivity.csv")
    reg, covs_used = regressions(dm, REG_VARS, covs, ref_ring=cfg.far_ring)
    numbers["covariates_used"] = covs_used
    _csv(reg, T / "Table8_regressions_cluster_robust.csv")
    df["coast_band"] = pd.cut(df["dist_coast_km"], [0, 10, 20, 50, 100, 150, np.inf], right=False).astype(str)
    strat = pd.DataFrame([stratified_contrast(df, v, cfg.near_ring, cfg.far_ring, sc, cfg.n_boot, cfg.random_seed)
                          for sc in (["coast_band"], ["coast_band", "UF"]) for v in CONTRAST_VARS])
    _csv(strat, T / "TableS18_coast_stratified_contrasts.csv")
    lobo = pd.concat([leave_one_block_out(df, v, cfg.near_ring, cfg.far_ring) for v in CONTRAST_VARS],
                     ignore_index=True)
    _csv(lobo, T / "TableS16_leave_one_block_out.csv")
    numbers["near_ring_effective_n_blocks"] = float(lobo["effective_n_near_blocks"].iloc[0]) if len(lobo) else None
    numbers["leave_one_block_out_range"] = {
        v: [float(g["estimate_without_block"].min()), float(g["estimate_without_block"].max())]
        for v, g in lobo.groupby("variable")}
    df.to_csv(out / "pixels_strict.csv.gz", index=False, compression="gzip")

    # block-size sensitivity (bootstrap and cluster-robust errors)
    brows, rrows = [], []
    for bdeg in cfg.sensitivity_blocks_deg:
        dfx = df.assign(block=block_ids(df["lon"], df["lat"], bdeg))
        for v in ("enh", "anom_sd", "tail_up_90", "tail_asym_90"):
            brows.append(block_bootstrap_compare(dfx, v, cfg.near_ring, cfg.far_ring, cfg.n_boot, cfg.random_seed)
                         .assign(block_deg=bdeg))
        dmx = dm.assign(block=block_ids(dm["lon"], dm["lat"], bdeg))
        rrows.append(regressions(dmx, ["enh", "anom_sd"], covs, ref_ring=cfg.far_ring)[0].assign(block_deg=bdeg))
    _csv(pd.concat(brows, ignore_index=True), T / "TableS10_block_size_bootstrap.csv")
    _csv(pd.concat(rrows, ignore_index=True), T / "TableS11_block_size_regressions.csv")
    # the v3 statistic, recomputed correctly (exclusive rings, same composite)
    sd_ring = results["relaxed"]["rings"].query("variable == 'xch4_mean'").set_index("ring")["sd"]
    numbers["v3_style_spatial_sd_by_ring_relaxed"] = {str(k): float(v) for k, v in sd_ring.items()}

    # ------------------------------------------------------------ time series
    bg_all = read_background_csv(gee_dir / "monthly_background.csv")
    trend, bg_fit = fit_trend(bg_all, cfg.min_bg_cells_per_month)
    numbers["regional_background_trend"] = trend
    _csv(bg_all.drop(columns="date"), T / "TableS5_monthly_background.csv")
    F.fig_background(bg_fit, trend, bg_all, cfg, FG)
    ts = None
    if len(ct) and has_monthly:
        MA, _ = read_multiband(gee_dir / "monthly_anom.tif")
        MC, _ = read_multiband(gee_dir / "monthly_count.tif")
        ma = {k[2:].replace("_", "-"): v for k, v in MA.items()}
        mc_ = {k[2:].replace("_", "-"): v for k, v in MC.items()}
        ts = cluster_monthly_series(labels, ma, mc_, ct["cluster"].head(6).tolist())
        _csv(ts.drop(columns="date"), T / "TableS6_cluster_monthly_series.csv")

    # ---------------------------------------------------------- sensitivity
    sens = None
    if sensitivity:
        rows, sig_rows = [], []
        mo_, mm_ = cfg.min_obs, cfg.min_valid_months   # rules scale with the run's own minimum
        count_rules = [(1, 0), (mo_, 0), (int(round(mo_ / 2)), mm_), (mo_, mm_), (2 * mo_, mm_), (4 * mo_, mm_)]
        for coastal in (True, False):
            st = strict_s if coastal else relaxed_s
            bgs = bg_s if coastal else has & (land >= 0.5)
            for mo, mm in count_rules:
                for zt in (2.0, 2.5, 3.0):
                    d = detect(f, grid, st, bgs, cfg, z_threshold=zt, min_obs=mo, min_months=mm)
                    pd_, _ = peak_distances(d, pipe_d)
                    scr = d["screen"]
                    rows.append(dict(coastal_screen=coastal, min_obs=mo, min_months=mm, z_threshold=zt,
                                     n_clusters=d["n_clusters"], n_negative_mirror=d["n_negative_mirror"],
                                     n_clusters_q=d["n_q"], n_within_5km=int(np.sum(pd_ <= 5)),
                                     n_within_20km=int(np.sum(pd_ <= 20)),
                                     # power context: how many screened cells are near a pipeline at all
                                     n_cells_screen=int(scr.sum()),
                                     n_cells_within_5km=int((scr & (pipe_d <= 5)).sum()),
                                     n_cells_within_20km=int((scr & (pipe_d <= 20)).sum())))
                    if d["n_q"]:
                        # chance that none of the n_q significant peaks falls within 90 km of a pipeline
                        # if peaks were placed at random over this rule's screened cells
                        p90 = float((scr & (pipe_d < 90)).sum() / max(1, scr.sum()))
                        rows[-1].update(share_cells_within_90km=p90,
                                        p_no_sig_peak_within_90km_random=float((1 - p90) ** d["n_q"]))
                        sig_rows.append(significant_clusters(d, grid, pipe_d, coast_d, f["n_obs"], cfg).assign(
                            coastal_screen=coastal, min_obs=mo, min_months=mm, z_threshold=zt))
        sens = pd.DataFrame(rows)
        _csv(sens, T / "TableS7_sensitivity.csv")
        _csv(pd.concat(sig_rows, ignore_index=True) if sig_rows else pd.DataFrame(columns=[
            "cluster", "peak_lon", "peak_lat", "peak_z", "q_mirror", "dist_peak_to_pipeline_km",
            "dist_peak_to_coast_km", "peak_n_obs", "coastal_screen", "min_obs", "min_months", "z_threshold"]),
             T / "TableS17_sensitivity_significant_clusters.csv")
        wrows = []
        for wk in cfg.sensitivity_windows_km:
            d = detect(f, grid, strict_s, bg_s, cfg, window_km=wk)
            wrows.append(dict(window_km=wk, n_clusters=d["n_clusters"], n_clusters_q=d["n_q"],
                              n_negative_mirror=d["n_negative_mirror"], robust_sigma_ppb=d["mad"]))
        _csv(pd.DataFrame(wrows), T / "TableS8_sensitivity_background_window.csv")

    # --------------------------------------------------------------- figures
    f_plot = dict(f, z=z, in_ne=in_ne)
    F.fig_study_area(ne, states, infra, lines, cfg, FG)
    F.fig_composite_and_density(f, grid, ne, lines, relaxed, cfg, FG)
    F.fig_anomaly_maps(f_plot, grid, ne, lines, polys, strict, cfg, FG, cfg.z_threshold)
    F.fig_histograms(legacy_x, mu if np.isfinite(mu) else np.nanmean(legacy_x),
                     sd if np.isfinite(sd) else np.nanstd(legacy_x),
                     enh[strict], D["med"], D["mad"], cfg.z_threshold, cfg, FG)
    F.fig_states(results["strict"]["df"], cfg, FG)
    F.fig_rings(results["strict"]["rings"], results["strict"]["order"], cfg, FG, panels=[
        ("anom_mean", "De-trended anomaly (ppb)"), ("enh", "Local enhancement (ppb)"),
        ("anom_sd", "Temporal SD (ppb)"), ("tail_up_90", "Upper tail, p90 − p50 (ppb)"),
        ("tail_lo_10", "Lower tail, p50 − p10 (ppb)"), ("dist_coast_km", "Distance to coast (km)")])
    if ts is not None:
        F.fig_cluster_timeseries(ts, ct["cluster"].tolist(), cfg, FG)
    if sens is not None:
        F.fig_sensitivity(sens, cfg, FG)

    # ---------------------------------------------------------------- report
    numbers["provenance"] = {k: provenance.get(k) for k in (
        "created_utc", "collection", "band", "start_date", "end_date", "run_label", "earthengine_api",
        "collection_versions", "run_properties", "assets", "files")}
    numbers["outputs"] = sorted(str(p.relative_to(out)) for p in out.rglob("*")
                                if p.is_file() and "gee" not in p.relative_to(out).parts)
    (out / "manuscript_numbers.json").write_text(json.dumps(_json_safe(numbers), indent=2, allow_nan=False),
                                                   encoding="utf-8")
    (out / "checks.json").write_text(json.dumps(_json_safe(checks), indent=2, allow_nan=False), encoding="utf-8")
    return dict(numbers=numbers, checks=checks, clusters=ct, enrichment=enrich, out=out)


def compare_runs(out_a: Path, out_b: Path, dest: Path):
    """Cells flagged as clusters in two runs (e.g. P1 and P2 processor periods).
    Writes the overlap per cluster of run A; a robust anomaly should persist."""
    A, ga = read_multiband(Path(out_a) / "derived_fields.tif")
    B, gb = read_multiband(Path(out_b) / "derived_fields.tif")
    if (ga.width, ga.height, tuple(ga.transform)) != (gb.width, gb.height, tuple(gb.transform)):
        raise RuntimeError("The two runs are not on the same grid")
    la = np.nan_to_num(A["cluster"]).astype(int)
    lb = np.nan_to_num(B["cluster"]).astype(int)
    rows = []
    for k in range(1, la.max() + 1):
        m = la == k
        rows.append(dict(cluster_run_a=k, n_pixels=int(m.sum()), n_pixels_also_cluster_in_b=int((m & (lb > 0)).sum()),
                         clusters_b_overlapping=",".join(str(x) for x in sorted(set(lb[m & (lb > 0)].tolist())))))
    df = pd.DataFrame(rows)
    df["persists"] = df["n_pixels_also_cluster_in_b"] > 0
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    _csv(df, dest)
    # near-minus-far contrasts on the cells that pass the strict screen in BOTH
    # runs, so a change between runs is not a change of sample composition
    pa = pd.read_csv(Path(out_a) / "pixels_strict.csv.gz")
    pb = pd.read_csv(Path(out_b) / "pixels_strict.csv.gz")
    key = ["lon", "lat"]
    common = pa[key].round(6).merge(pb[key].round(6), on=key)
    crow = []
    for tag, px in (("a", pa), ("b", pb)):
        px = px.assign(lon=px["lon"].round(6), lat=px["lat"].round(6)).merge(common, on=key)
        near, far = px["ring"] == "0-5 km", px["ring"] == ">=50 km"
        for v in ("anom_mean", "enh", "anom_sd", "tail_up_90", "tail_lo_10", "tail_asym_90"):
            crow.append(dict(run=tag, variable=v, near_minus_far=float(px.loc[near, v].mean() - px.loc[far, v].mean()),
                             n_near_common=int(near.sum()), n_far_common=int(far.sum())))
    _csv(pd.DataFrame(crow), Path(dest).with_name(Path(dest).stem.replace("TableS12_clusters", "TableS19_near_far_common_cells") + ".csv"))
    return df
