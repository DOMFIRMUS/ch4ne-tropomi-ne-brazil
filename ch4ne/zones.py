"""
Proximity-to-infrastructure analysis.

The v3 "variance signature" compared the SPATIAL standard deviation of a
multi-year mean inside nested buffers (5 km within 15 km within the whole
region) and read it as TEMPORAL intermittency. Pipelines in NE Brazil run
along the coast, where land/sea-mixed cells and dark surfaces widen the
distribution in BOTH tails. This module therefore:

1. uses mutually exclusive distance rings;
2. reports, per ring, the temporal variability (anom_sd), the widths of the
   upper AND lower tails of each cell's anomaly distribution (tail_up /
   tail_lo; intermittent emissions should widen only the upper tail), and
   robust spread measures with the upper and lower semi-deviations
   separately;
3. compares near-pipeline pixels with covariate-matched far pixels (distance
   to coast, latitude, albedo, water, elevation, data density);
4. fits regressions with cluster-robust (spatial block) standard errors,
   including a dispersion model on absolute residuals;
5. quantifies uncertainty by a spatial block bootstrap.
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from scipy.spatial import cKDTree

COVARIATES = ["dist_coast_km", "lat", "albedo_b7", "water_frac", "elev_mean", "elev_sd",
              "log_n_obs", "frac_janjun", "aod_mean"]


def ring_labels(edges):
    labs = []
    for a, b in zip(edges[:-1], edges[1:]):
        labs.append(f">={a:g} km" if np.isinf(b) else f"{a:g}-{b:g} km")
    return labs


def assign_rings(dist_km, edges):
    labs = ring_labels(edges)
    out = np.full(np.shape(dist_km), None, dtype=object)
    for lab, a, b in zip(labs, edges[:-1], edges[1:]):
        out[(dist_km >= a) & (dist_km < b)] = lab
    return out, labs


def semi_sd(x):
    """Upper and lower semi-deviations about the median."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    m = np.median(x)
    up = x[x > m] - m
    lo = m - x[x < m]
    return (np.sqrt(np.mean(up ** 2)) if up.size else np.nan,
            np.sqrt(np.mean(lo ** 2)) if lo.size else np.nan)


def describe(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return dict(n=int(x.size))
    q = np.percentile(x, [5, 25, 50, 75, 95])
    su, sl = semi_sd(x)
    return dict(n=int(x.size), mean=x.mean(), sd=x.std(ddof=1), median=q[2], q25=q[1], q75=q[3],
                iqr=q[3] - q[1], mad_sigma=1.4826 * np.median(np.abs(x - q[2])),
                semi_sd_upper=su, semi_sd_lower=sl, skew=stats.skew(x), p05=q[0], p95=q[4])


def ring_table(df, variables, ring_col="ring", order=None):
    rows = []
    for v in variables:
        for ring, g in df.groupby(ring_col, sort=False):
            d = describe(g[v])
            d.update(variable=v, ring=ring)
            rows.append(d)
    out = pd.DataFrame(rows)
    if order is not None:
        out["ring"] = pd.Categorical(out["ring"], order, ordered=True)
        out = out.sort_values(["variable", "ring"])
    cols = ["variable", "ring"] + [c for c in out.columns if c not in ("variable", "ring")]
    return out[cols].reset_index(drop=True)


def block_ids(lon, lat, block_deg):
    bx = np.floor(np.asarray(lon) / block_deg).astype(int)
    by = np.floor(np.asarray(lat) / block_deg).astype(int)
    return bx * 10000 + by


def _stat_pair(a, b):
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return dict(mean_diff=np.nan, median_diff=np.nan, sd_ratio=np.nan, iqr_ratio=np.nan,
                    upper_semi_sd_ratio=np.nan, lower_semi_sd_ratio=np.nan)
    ua, la = semi_sd(a)
    ub, lb = semi_sd(b)
    iqr = lambda x: np.subtract(*np.percentile(x, [75, 25]))
    return dict(mean_diff=a.mean() - b.mean(), median_diff=np.median(a) - np.median(b),
                sd_ratio=a.std(ddof=1) / b.std(ddof=1), iqr_ratio=iqr(a) / iqr(b),
                upper_semi_sd_ratio=ua / ub, lower_semi_sd_ratio=la / lb)


def block_bootstrap_compare(df, var, near, far, n_boot, seed, ring_col="ring", block_col="block"):
    """Near-vs-far contrasts with 95 % CIs from resampling spatial blocks.

    Blocks (default 1 x 1 deg) are drawn with replacement; all pixels of a drawn
    block enter with their own ring label, which preserves within-block spatial
    correlation (Kunsch 1989; Lahiri 2003)."""
    sub = df[df[ring_col].isin([near, far])][[var, ring_col, block_col]].dropna()
    point = _stat_pair(sub.loc[sub[ring_col] == near, var].values, sub.loc[sub[ring_col] == far, var].values)
    blocks = sub[block_col].unique()
    groups = {b: g for b, g in sub.groupby(block_col)}
    rng = np.random.default_rng(seed)
    reps = []
    for _ in range(n_boot if blocks.size >= 5 else 0):   # no interval from fewer than 5 blocks
        pick = rng.choice(blocks, size=blocks.size, replace=True)
        s = pd.concat([groups[b] for b in pick], ignore_index=True)
        reps.append(_stat_pair(s.loc[s[ring_col] == near, var].values, s.loc[s[ring_col] == far, var].values))
    reps = pd.DataFrame(reps, columns=list(point))
    rows = []
    for k, v in point.items():
        lo, hi = np.nanpercentile(reps[k], [2.5, 97.5]) if len(reps) else (np.nan, np.nan)
        rows.append(dict(variable=var, contrast=k, near=near, far=far, estimate=v, ci95_low=lo, ci95_high=hi,
                         n_near=int((sub[ring_col] == near).sum()), n_far=int((sub[ring_col] == far).sum()),
                         n_blocks=int(blocks.size)))
    bf = stats.levene(sub.loc[sub[ring_col] == near, var], sub.loc[sub[ring_col] == far, var], center="median")
    rows.append(dict(variable=var, contrast="brown_forsythe_p_naive", near=near, far=far,
                     estimate=bf.pvalue, ci95_low=np.nan, ci95_high=np.nan,
                     n_near=rows[0]["n_near"], n_far=rows[0]["n_far"], n_blocks=int(blocks.size)))
    return pd.DataFrame(rows)


def propensity_logit(df, near_mask, far_mask, covariates):
    """Linear predictor (logit) of the probability of being a near-pipeline
    pixel, from a logistic regression on the standardised covariates fitted to
    the near and far pixels. NaN outside near|far. A small ridge penalty is
    used only if the unpenalised fit does not converge (quasi-separation)."""
    import statsmodels.api as sm
    sel = np.asarray(near_mask | far_mask)
    X = df.loc[sel, covariates]
    Z = ((X - X.mean()) / X.std().replace(0, 1)).values
    A = sm.add_constant(Z, has_constant="add")
    y = np.asarray(near_mask)[sel].astype(float)
    method = "logit"
    try:
        fit = sm.Logit(y, A).fit(disp=0, maxiter=500)
        if not fit.mle_retvals.get("converged", True) or not np.all(np.isfinite(fit.params)):
            raise RuntimeError("not converged")
        params = fit.params
    except Exception:
        params = sm.GLM(y, A, family=sm.families.Binomial()).fit_regularized(alpha=1e-3, L1_wt=0.0).params
        method = "logit_ridge"
    out = np.full(len(df), np.nan)
    out[sel] = A @ params
    return out, method


def matched_controls(df, near_mask, far_mask, covariates, k, caliper, balance_extra=(),
                     method="ps_mahalanobis"):
    """For every near-pipeline pixel, the k nearest far pixels (with
    replacement), keeping only matches within the caliper.

    method="ps_mahalanobis" (default): Mahalanobis distance on the covariates
        among the far pixels whose propensity-score logit lies within the
        caliper (in SDs of the logit; 0.2 as recommended by Austin 2011,
        Pharm. Stat. 10:150) - Rosenbaum & Rubin (1985, Am. Stat. 39:33).
    method="propensity": nearest on the propensity-score logit alone.
    method="euclidean": distance in the standardised covariate space, caliper
        in standardised units (strict in many dimensions; sensitivity only).

    Returns (controls, balance table, mask of matched near pixels, summary).
    `controls` rows are copies of the far pixels with a 'match_of' column.
    `balance_extra` lists variables reported in the balance table without
    being matched on (e.g. the sampling time mean_t_years)."""
    near_mask, far_mask = np.asarray(near_mask), np.asarray(far_mask)
    ps_fit = None
    near_idx = np.where(near_mask)[0]
    far_idx = np.where(far_mask)[0]
    if method == "ps_mahalanobis":
        lp, ps_fit = propensity_logit(df, near_mask, far_mask, covariates)
        caliper_abs = caliper * _pooled_sd(lp, near_mask, far_mask)
        X = df[covariates].values
        pool = near_mask | far_mask
        VI = np.linalg.pinv(np.cov(X[pool], rowvar=False))
        order = np.argsort(lp[far_idx])
        lp_far = lp[far_idx][order]
        dist = np.full((near_idx.size, k), np.inf)
        j = np.zeros((near_idx.size, k), dtype=int)
        for r, i in enumerate(near_idx):
            lo = np.searchsorted(lp_far, lp[i] - caliper_abs, side="left")
            hi = np.searchsorted(lp_far, lp[i] + caliper_abs, side="right")
            cand = order[lo:hi]
            if cand.size == 0:
                continue
            d = X[far_idx[cand]] - X[i]
            md = np.sqrt(np.einsum("ij,jk,ik->i", d, VI, d))
            take = np.argsort(md)[:k]
            dist[r, :take.size] = md[take]
            j[r, :take.size] = cand[take]
        ok = np.isfinite(dist)
        return _finish_matching(df, near_mask, far_mask, covariates, k, caliper, caliper_abs, method, ps_fit,
                                near_idx, far_idx, dist, j, ok, balance_extra, lp)
    if method == "propensity":
        lp, ps_fit = propensity_logit(df, near_mask, far_mask, covariates)
        Z = lp[:, None]
        caliper_abs = caliper * _pooled_sd(lp, near_mask, far_mask)
    elif method == "euclidean":
        X = df[covariates]
        mu, sd = X[far_mask | near_mask].mean(), X[far_mask | near_mask].std().replace(0, 1)
        Z = ((X - mu) / sd).values
        caliper_abs = caliper
    else:
        raise ValueError(method)
    tree = cKDTree(Z[far_idx])
    dist, j = tree.query(Z[near_idx], k=k)
    dist, j = np.atleast_2d(dist), np.atleast_2d(j)
    if j.shape[0] != near_idx.size:
        dist, j = dist.T, j.T
    ok = dist <= caliper_abs
    return _finish_matching(df, near_mask, far_mask, covariates, k, caliper, caliper_abs, method, ps_fit,
                            near_idx, far_idx, dist, j, ok, balance_extra, lp if method == "propensity" else None)


def _pooled_sd(x, a_mask, b_mask):
    """sqrt((s_a^2 + s_b^2) / 2), the pooled within-group SD of Austin (2011)."""
    return float(np.sqrt((np.nanvar(x[a_mask], ddof=1) + np.nanvar(x[b_mask], ddof=1)) / 2))


def _finish_matching(df, near_mask, far_mask, covariates, k, caliper, caliper_abs, method, ps_fit,
                     near_idx, far_idx, dist, j, ok, balance_extra, lp=None):
    rows_near = np.repeat(near_idx, k).reshape(j.shape)
    ctrl = df.iloc[far_idx[j[ok]]].copy()
    ctrl["match_of"] = df.index[rows_near[ok]].values
    ctrl["match_dist"] = dist[ok]
    matched = np.zeros(len(df), bool)
    matched[near_idx[ok.any(axis=1)]] = True
    # each matched near pixel counts once; its m_i controls share its weight
    # (weight 1/m_i per control row), as in the matching estimator itself
    ctrl["w"] = 1.0 / ctrl.groupby("match_of")["match_of"].transform("size")
    cols = list(covariates) + [c for c in balance_extra if c not in covariates]
    vals = {c: df[c].values.astype(float) for c in cols}
    if lp is not None:
        cols.append("ps_logit")
        vals["ps_logit"] = lp
    bal = []
    for c in cols:
        v = vals[c]
        near_all, near_m, far_all = v[near_mask], v[matched], v[far_mask]
        # fixed denominator: pooled SD of ALL near and ALL far cells (before matching)
        pooled = np.sqrt((np.nanvar(near_all, ddof=1) + np.nanvar(far_all, ddof=1)) / 2)
        sd_near = np.nanstd(near_all, ddof=1)
        pooled = pooled if np.isfinite(pooled) and pooled > 0 else np.nan
        mm = float(np.average(v[df.index.get_indexer(ctrl.index)], weights=ctrl["w"])) if len(ctrl) else np.nan
        mn = float(np.nanmean(near_m)) if near_m.size else np.nan
        bal.append(dict(covariate=c, matched_on=c in covariates, mean_near=float(np.nanmean(near_all)),
                        mean_near_matched=mn, mean_far_all=float(np.nanmean(far_all)), mean_matched_controls=mm,
                        smd_before=(np.nanmean(near_all) - np.nanmean(far_all)) / pooled,
                        smd_after=(mn - mm) / pooled,
                        smd_after_near_sd=(mn - mm) / sd_near if np.isfinite(sd_near) and sd_near > 0 else np.nan))
    bal = pd.DataFrame(bal)
    used = ctrl.index.value_counts()
    wsum = ctrl.groupby(level=0)["w"].sum().sort_values(ascending=False) if len(ctrl) else pd.Series(dtype=float)
    on = bal[bal.matched_on]
    smd_ok = len(on) and on["smd_after"].notna().all()
    sizes = ctrl.groupby("match_of").size().value_counts().sort_index() if len(ctrl) else pd.Series(dtype=int)
    summary = dict(method=method, propensity_fit=ps_fit, caliper=caliper, caliper_abs=float(caliper_abs), k=k,
                   n_near=int(near_mask.sum()), n_near_matched=int(matched.sum()),
                   n_near_unmatched=int(near_mask.sum() - matched.sum()), n_control_rows=int(len(ctrl)),
                   n_distinct_controls=int(used.size), max_reuse_of_one_control=int(used.max()) if len(used) else 0,
                   # effective number of distinct controls (Kish) given the 1/m_i weights and reuse
                   effective_n_controls=float(wsum.sum() ** 2 / (wsum ** 2).sum()) if len(wsum) else 0.0,
                   weight_share_top3_controls=float(wsum.head(3).sum() / wsum.sum()) if len(wsum) else float("nan"),
                   top3_controls=[dict(lon=float(df.loc[i, "lon"]), lat=float(df.loc[i, "lat"]),
                                       UF=str(df.loc[i, "UF"]) if "UF" in df else None,
                                       dist_coast_km=float(df.loc[i, "dist_coast_km"]) if "dist_coast_km" in df else None,
                                       weight_share=float(wsum.loc[i] / wsum.sum()))
                                  for i in wsum.head(3).index] if len(wsum) else [],
                   control_set_sizes={int(kk): int(vv) for kk, vv in sizes.items()},
                   max_abs_smd_after=float(on["smd_after"].abs().max()) if smd_ok else float("nan"),
                   n_covariates_abs_smd_after_gt_0_1=int((on["smd_after"].abs() > 0.1).sum()) if smd_ok else None)
    if lp is not None:  # common support of the propensity score
        far_hi = np.nanpercentile(lp[far_mask], 99)
        summary.update(ps_logit_median_near=float(np.nanmedian(lp[near_mask])),
                       ps_logit_p99_far=float(far_hi),
                       share_near_above_far_p99=float(np.mean(lp[near_mask] > far_hi)))
    return ctrl, bal, matched, summary


def _wquantile(x, w, q):
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(q, c, x)


def _wsemi_sd(x, w):
    m = _wquantile(x, w, 0.5)
    up, lo = x > m, x < m
    return (np.sqrt(np.average((x[up] - m) ** 2, weights=w[up])) if up.any() else np.nan,
            np.sqrt(np.average((m - x[lo]) ** 2, weights=w[lo])) if lo.any() else np.nan)


def _wmedian(x, w):
    return _wquantile(x, w, 0.5)


def _wsd(x, w, freq):
    m = np.average(x, weights=w)
    v = np.average((x - m) ** 2, weights=w)
    if freq:   # frequency weights (bootstrap multiplicities)
        return np.sqrt(v * w.sum() / (w.sum() - 1))
    return np.sqrt(v * w.sum() ** 2 / (w.sum() ** 2 - (w ** 2).sum()))   # reliability weights


def _matched_stats(y_near, w_near, set_means: dict, y_ctrl, w_ctrl):
    """Matching estimator: (weighted) mean/median over matched near pixels of
    Y_i minus the mean of ITS OWN control set (Abadie & Imbens 2006,
    Econometrica 74:235), for the raw set means and for each bias-corrected
    version in `set_means`; distribution ratios compare the near values with
    the control values weighted 1/m_i, so every near pixel's control set
    carries equal weight. w_near are multiplicities (1 in the point estimate)."""
    out = dict(mean_diff=np.nan, median_diff=np.nan, sd_ratio=np.nan, iqr_ratio=np.nan,
               upper_semi_sd_ratio=np.nan, lower_semi_sd_ratio=np.nan)
    for name, sm in set_means.items():
        ok = np.isfinite(y_near) & np.isfinite(sm) & (w_near > 0)
        suffix = "" if name == "raw" else f"_{name}"
        if ok.sum() >= 3:
            d, w = y_near[ok] - sm[ok], w_near[ok]
            out[f"mean_diff{suffix}"] = float(np.average(d, weights=w))
            out[f"median_diff{suffix}"] = float(_wmedian(d, w))
        else:
            out[f"mean_diff{suffix}"] = out[f"median_diff{suffix}"] = np.nan
    ga = np.isfinite(y_near) & (w_near > 0)
    gb = np.isfinite(y_ctrl) & (w_ctrl > 0)
    a, wa, b, wb = y_near[ga], w_near[ga], y_ctrl[gb], w_ctrl[gb]
    if wa.sum() >= 3 and b.size >= 3:
        iqr_a = _wquantile(a, wa, 0.75) - _wquantile(a, wa, 0.25)
        iqr_b = _wquantile(b, wb, 0.75) - _wquantile(b, wb, 0.25)
        ua, la = _wsemi_sd(a, wa)
        ub, lb = _wsemi_sd(b, wb)
        sd_b = _wsd(b, wb, False)
        div = lambda x, y: x / y if np.isfinite(y) and y > 0 else np.nan
        out.update(sd_ratio=div(_wsd(a, wa, True), sd_b), iqr_ratio=div(iqr_a, iqr_b),
                   upper_semi_sd_ratio=div(ua, ub), lower_semi_sd_ratio=div(la, lb))
    return out


def _wls_beta(X, y, w):
    A = np.column_stack([np.ones(len(y)), X]) * np.sqrt(w)[:, None]
    return np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0][1:]


def matched_contrast(df, ctrl, matched_mask, var, n_boot, seed, block_col="block", covariates=None, far_mask=None,
                     min_blocks=5):
    """Matched near pixels vs their own control sets.

    Point estimate: the matching estimator of the effect on the near pixels,
    mean_i [Y_i - mean_{j in J(i)} Y_j] (each near pixel counts once, however
    many controls it kept). With `covariates`, also bias-corrected versions
    (Abadie & Imbens 2011, J. Bus. Econ. Stat. 29:1): each control's value is
    shifted by beta'(x_i - x_j), with beta from a regression of `var` on the
    covariates over the MATCHED controls weighted by their use (primary,
    '_bias_corrected'; it avoids extrapolating from far cells without common
    support) and, as a sensitivity, over the whole far pool
    ('_bias_corrected_farpool').
    Uncertainty: joint spatial block bootstrap - blocks are drawn from the
    union of the near and control blocks; near pixels and controls enter with
    their own block's multiplicity and every set mean is recomputed from the
    drawn controls (sets left without controls drop out); matches and beta
    are held fixed. With fewer than `min_blocks` blocks no interval is given."""
    names = ["mean_diff", "median_diff", "sd_ratio", "iqr_ratio", "upper_semi_sd_ratio", "lower_semi_sd_ratio"]
    has_adj = covariates is not None
    if has_adj:
        names += ["mean_diff_bias_corrected", "median_diff_bias_corrected"]
        if far_mask is not None:
            names += ["mean_diff_bias_corrected_farpool", "median_diff_bias_corrected_farpool"]
    near_ids = df.index[np.asarray(matched_mask)]
    if len(near_ids) == 0 or len(ctrl) == 0:  # nothing matched (e.g. strict caliper, short period)
        return pd.DataFrame([dict(variable=var, contrast=k, estimate=np.nan, ci95_low=np.nan, ci95_high=np.nan,
                                  n_near=0, n_control_rows=int(len(ctrl)), n_blocks=0) for k in names])
    y_near = df.loc[near_ids, var].values.astype(float)
    bn = df.loc[near_ids, block_col].values
    bc = ctrl[block_col].values
    y_c = ctrl[var].values.astype(float)
    w_row = ctrl["w"].values if "w" in ctrl else \
        (1.0 / ctrl.groupby("match_of")["match_of"].transform("size")).values
    pos_of = {nid: i for i, nid in enumerate(near_ids)}
    c_pos = ctrl["match_of"].map(pos_of).values.astype(int)
    nN = len(near_ids)
    adj = {"raw": y_c}
    if has_adj:
        xn = df.loc[ctrl["match_of"].values, covariates].values
        xc = ctrl[covariates].values
        wj = ctrl.groupby(level=0)["w"].sum() if "w" in ctrl else ctrl.groupby(level=0).size().astype(float)
        dc = df.loc[wj.index, [var] + list(covariates)].dropna()
        A_ = np.column_stack([np.ones(len(dc)), dc[covariates].values])
        # enough distinct controls (>= 5 per parameter) and full rank, else no correction
        if len(dc) >= 5 * A_.shape[1] and np.linalg.matrix_rank(A_) == A_.shape[1]:
            beta_m = _wls_beta(dc[covariates].values, dc[var].values, wj.loc[dc.index].values)
            adj["bias_corrected"] = y_c + (xn - xc) @ beta_m
        else:
            adj["bias_corrected"] = np.full(len(y_c), np.nan)
        if far_mask is not None:
            pool = df.loc[far_mask, [var] + list(covariates)].dropna()
            beta_f = _wls_beta(pool[covariates].values, pool[var].values, np.ones(len(pool)))
            adj["bias_corrected_farpool"] = y_c + (xn - xc) @ beta_f

    def set_means(wc):
        sw = np.bincount(c_pos, weights=wc, minlength=nN)
        out = {}
        for k2, yy in adj.items():
            g = np.isfinite(yy)
            sy = np.bincount(c_pos[g], weights=(wc * yy)[g], minlength=nN)
            swg = np.bincount(c_pos[g], weights=wc[g], minlength=nN)
            with np.errstate(invalid="ignore", divide="ignore"):
                out[k2] = np.where(swg > 0, sy / np.where(swg > 0, swg, 1), np.nan)
        return out, sw

    ones_c = np.ones(len(y_c))
    sm0, _ = set_means(ones_c)
    point = _matched_stats(y_near, np.ones(nN), sm0, y_c, w_row)
    blocks = np.unique(np.concatenate([bn, bc]))
    bpos = {b_: i for i, b_ in enumerate(blocks)}
    ibn = np.array([bpos[b_] for b_ in bn], dtype=int)
    ibc = np.array([bpos[b_] for b_ in bc], dtype=int)
    rows_ci = {k2: (np.nan, np.nan) for k2 in point}
    valid = {k2: 0 for k2 in point}
    n_near_blocks = int(np.unique(bn).size)
    if n_near_blocks >= min_blocks:          # the near cells carry the information
        rng = np.random.default_rng(seed)
        reps = []
        for _ in range(n_boot):
            cnt = np.bincount(rng.integers(0, blocks.size, blocks.size), minlength=blocks.size).astype(float)
            wn, wc = cnt[ibn], cnt[ibc]
            sm, sw = set_means(wc)
            wn_eff = np.where(sw > 0, wn, 0.0)
            # each drawn set carries its near pixel's multiplicity, split over its drawn controls
            w_c = wc / np.where(sw[c_pos] > 0, sw[c_pos], 1.0) * wn_eff[c_pos]
            reps.append(_matched_stats(y_near, wn_eff, sm, y_c, w_c))
        reps = pd.DataFrame(reps)
        for k2 in point:
            ok = np.isfinite(reps[k2].values)
            valid[k2] = int(ok.sum())
            # no interval if more than 5 % of the replicates are undefined
            if ok.mean() >= 0.95:
                rows_ci[k2] = tuple(np.percentile(reps[k2].values[ok], [2.5, 97.5]))
    rows = []
    for k2, v in point.items():
        lo, hi = rows_ci[k2]
        rows.append(dict(variable=var, contrast=k2, estimate=v, ci95_low=lo, ci95_high=hi,
                         n_near=int(nN), n_control_rows=int(len(ctrl)), n_blocks=int(blocks.size),
                         n_near_blocks=n_near_blocks, n_valid_reps=valid[k2]))
    return pd.DataFrame(rows)


def stratified_contrast(df, var, near, far, strata_cols, n_boot, seed, ring_col="ring", block_col="block",
                        min_blocks=5):
    """Near-minus-far mean difference within exact strata (e.g. coast-distance
    band, or band x state), averaged with the near-cell counts as weights
    (the effect on the near cells). Near cells in strata without far cells are
    left out and counted. Joint spatial block bootstrap for the CI."""
    cols = [var, ring_col, block_col] + list(strata_cols)
    sub = df[df[ring_col].isin([near, far])][cols].dropna()
    key = sub[list(strata_cols)].astype(str).agg("|".join, axis=1).values
    is_n = (sub[ring_col] == near).values
    y = sub[var].values.astype(float)
    blk = sub[block_col].values

    def est(w):
        tot, acc, used = 0.0, 0.0, 0.0
        for k2 in np.unique(key):
            m = key == k2
            wn, wf = w[m & is_n], w[m & ~is_n]
            if wn.sum() > 0 and wf.sum() > 0:
                d = np.average(y[m & is_n], weights=wn) - np.average(y[m & ~is_n], weights=wf)
                acc += wn.sum() * d
                used += wn.sum()
            tot += wn.sum()
        return (acc / used if used else np.nan), used, tot

    point, used, tot = est(np.ones(len(y)))
    # strata that carry information: both near and far cells present
    informative = np.zeros(len(y), bool)
    for k2 in np.unique(key):
        m = key == k2
        if (m & is_n).any() and (m & ~is_n).any():
            informative |= m
    nb_near = int(np.unique(blk[informative & is_n]).size)
    nb_far = int(np.unique(blk[informative & ~is_n]).size)
    blocks = np.unique(blk)
    lo = hi = np.nan
    n_valid = 0
    if nb_near >= min_blocks:
        bpos = {b_: i for i, b_ in enumerate(blocks)}
        ib = np.array([bpos[b_] for b_ in blk], dtype=int)
        rng = np.random.default_rng(seed)
        reps = np.array([est(np.bincount(rng.integers(0, blocks.size, blocks.size),
                                         minlength=blocks.size)[ib].astype(float))[0] for _ in range(n_boot)])
        ok = np.isfinite(reps)
        n_valid = int(ok.sum())
        if ok.mean() >= 0.95:
            lo, hi = np.percentile(reps[ok], [2.5, 97.5])
    return dict(variable=var, strata="x".join(strata_cols), estimate=point, ci95_low=lo, ci95_high=hi,
                n_near=int(tot), n_near_in_strata_with_far=int(used), n_near_blocks_used=nb_near,
                n_far_blocks_used=nb_far, n_valid_reps=n_valid)


def leave_one_block_out(df, var, near, far, ring_col="ring", block_col="block"):
    """Near-minus-far mean difference recomputed with each spatial block that
    holds near pixels removed (both its near and far pixels), plus the
    effective number of near blocks (1 / Herfindahl index of the near-pixel
    shares). Shows how much a contrast rests on one or two blocks."""
    sub = df[df[ring_col].isin([near, far])][[var, ring_col, block_col]].dropna()
    is_n = sub[ring_col] == near
    full = sub.loc[is_n, var].mean() - sub.loc[~is_n, var].mean()
    share = sub.loc[is_n, block_col].value_counts(normalize=True)
    rows = []
    for blk, sh in share.items():
        s2 = sub[sub[block_col] != blk]
        n2 = s2[ring_col] == near
        est = s2.loc[n2, var].mean() - s2.loc[~n2, var].mean() if n2.any() else np.nan
        rows.append(dict(variable=var, dropped_block=int(blk), near_share_of_block=float(sh),
                         estimate_without_block=float(est), estimate_all=float(full)))
    out = pd.DataFrame(rows).sort_values("near_share_of_block", ascending=False)
    out["effective_n_near_blocks"] = float(1.0 / np.sum(share.values ** 2))
    return out


def regressions(df, y_vars, covariates, ring_col="ring", ref_ring=None, block_col="block"):
    """OLS of each response on ring dummies + covariates with cluster-robust SEs
    (clusters = spatial blocks), and a dispersion model.

    Dispersion model (Brown-Forsythe type): residuals of the FULL mean model
    (ring + covariates), centred on their ring median, in absolute value,
    regressed on ring dummies + covariates. Taking residuals from a model
    without the ring term would let a pure ring mean shift appear as extra
    dispersion."""
    rows = []
    d0 = df.dropna(subset=list(covariates) + [ring_col, block_col]).copy()
    covariates = drop_collinear(d0, covariates)
    ref = ref_ring or d0[ring_col].iloc[0]
    rhs_cov = " + ".join(covariates)
    ring_term = f'C({ring_col}, Treatment(reference="{ref}"))'
    for y in y_vars:
        d = d0.dropna(subset=[y])
        groups = pd.factorize(d[block_col])[0]  # statsmodels needs codes >= 0
        full = None
        for model, formula in (("mean", f"{y} ~ {ring_term} + {rhs_cov}"), ("mean_unadjusted", f"{y} ~ {ring_term}")):
            fit = smf.ols(formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": groups})
            rows += _ring_coefs(fit, y, model, ring_col)
            if model == "mean":
                full = fit
        r = full.resid
        d = d.assign(abs_dev=np.abs(r - r.groupby(d[ring_col]).transform("median")))
        fit = smf.ols(f"abs_dev ~ {ring_term} + {rhs_cov}", data=d).fit(
            cov_type="cluster", cov_kwds={"groups": groups})
        rows += _ring_coefs(fit, y, "dispersion_brown_forsythe", ring_col)
    return pd.DataFrame(rows), covariates


def drop_collinear(df, covariates, max_abs_corr=0.95):
    """Keep covariates in order, skipping any that is (nearly) constant or
    correlated above max_abs_corr with one already kept."""
    kept = []
    for c in covariates:
        x = df[c]
        if x.std() == 0:
            continue
        if all(abs(np.corrcoef(x, df[k])[0, 1]) <= max_abs_corr for k in kept):
            kept.append(c)
    return kept


def _ring_coefs(fit, y, model, ring_col):
    out = []
    ci = fit.conf_int()
    for name in fit.params.index:
        if name.startswith(f"C({ring_col}"):
            level = name.split("[T.")[-1].rstrip("]")
            out.append(dict(response=y, model=model, ring=level, coef=fit.params[name], se=fit.bse[name],
                            ci95_low=ci.loc[name, 0], ci95_high=ci.loc[name, 1], p=fit.pvalues[name],
                            n=int(fit.nobs), r2=fit.rsquared))
    return out


def coast_collocation(pipe_dist_km, coast_dist_km, region, near_km=5.0, coast_bands=(10.0, 20.0, 50.0)):
    """How strongly pipeline proximity coincides with coastal location.

    Uses geometry only (ANP lines, coastline, NE grid cells), so it can be
    computed before any satellite data are read."""
    near = region & (pipe_dist_km <= near_km)
    rows = []
    for name, m in (("cells <= %g km from a pipeline" % near_km, near), ("all NE cells", region)):
        d = coast_dist_km[m]
        row = dict(group=name, n_cells=int(m.sum()), median_dist_coast_km=float(np.median(d)))
        for b in coast_bands:
            row[f"share_within_{b:g}km_of_coast"] = float(np.mean(d <= b))
        rows.append(row)
    return pd.DataFrame(rows)

