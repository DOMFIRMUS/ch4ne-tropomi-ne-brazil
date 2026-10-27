"""
Temporal analyses: regional background growth and cluster time series.

v3 plotted RAW monthly medians at three polygon centroids, from a collection
that skipped the uncertainty filter, and read the 2020-2024 rise (~7.6 ppb/yr)
as local behaviour. Here the NE-wide monthly background is fitted with a
linear trend + annual and semi-annual harmonics (its slope is the regional
growth rate that must be removed), and cluster series are shown as
DE-TRENDED anomalies averaged over all cluster pixels, weighted by the number
of retrievals.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm


def read_background_csv(path):
    bg = pd.read_csv(path)
    cols = [c for c in ("month", "bg_ppb", "n_cells", "n_orbits", "lat_median_observed") if c in bg.columns]
    bg = bg[cols].copy()
    bg["date"] = pd.to_datetime(bg["month"] + "-15")
    bg = bg.sort_values("date").reset_index(drop=True)
    bg.loc[bg["bg_ppb"] <= -9990, "bg_ppb"] = np.nan
    return bg


def _design(t):
    return np.column_stack([np.ones(t.size), t,
                            np.sin(2 * np.pi * t), np.cos(2 * np.pi * t),
                            np.sin(4 * np.pi * t), np.cos(4 * np.pi * t)])


def _gap_aware_newey_west(X, resid, month_index, max_lag):
    """Newey-West covariance with Bartlett weights on the CALENDAR lag in
    months, so rows on either side of dropped months are not treated as
    adjacent (statsmodels' HAC counts lags in row order)."""
    n = X.shape[0]
    xe = X * resid[:, None]
    S = xe.T @ xe
    for i in range(n):
        for j in range(i + 1, n):
            lag = abs(month_index[j] - month_index[i])
            if lag <= max_lag:
                w = 1.0 - lag / (max_lag + 1.0)
                S += w * (np.outer(xe[i], xe[j]) + np.outer(xe[j], xe[i]))
    bread = np.linalg.inv(X.T @ X)
    return bread @ S @ bread


PROCESSOR_STEPS = ("2021-07-01", "2022-07-17")   # bias-correction changes in the OFFL product


def fit_trend(bg: pd.DataFrame, min_cells: int, lags=(3, 6), min_side_months=6):
    """bg_ppb ~ a + b*t + annual + semi-annual harmonics (OLS).

    Reports the slope with gap-aware Newey-West SEs for several lags, the
    peak-to-trough half-amplitude of the whole fitted seasonal cycle, and the
    number of months dropped for insufficient coverage."""
    ok = (bg["n_cells"] >= min_cells) & bg["bg_ppb"].notna()
    d = bg[ok].copy()
    t = ((d["date"] - pd.Timestamp("2020-01-01")).dt.days / 365.25).values
    X = _design(t)
    y = d["bg_ppb"].values
    fit = sm.OLS(y, X).fit()
    mi = (d["date"].dt.year * 12 + d["date"].dt.month).values
    se = {f"slope_se_nw_lag{L}": float(np.sqrt(_gap_aware_newey_west(X, fit.resid, mi, L)[1, 1])) for L in lags}
    ph = np.linspace(0, 1, 1000, endpoint=False)
    S = np.column_stack([np.sin(2 * np.pi * ph), np.cos(2 * np.pi * ph),
                         np.sin(4 * np.pi * ph), np.cos(4 * np.pi * ph)]) @ fit.params[2:6]
    d["fitted"] = fit.fittedvalues
    # same model plus level steps at the processor changes, to show how much
    # of the apparent growth may be a product-version artefact. The step dates
    # are imposed a priori (2 of the transitions between processor versions)
    # and are collinear with the trend, so the step model can absorb part of
    # the real growth; its slope and steps are reported with the same
    # gap-aware Newey-West errors. The month-by-month background subtraction
    # removes only the spatially uniform part of such steps.
    steps = np.column_stack([(d["date"] >= pd.Timestamp(s)).astype(float).values for s in PROCESSOR_STEPS])
    # a step is estimable only with enough months on BOTH sides; otherwise it
    # would just absorb the first or last month(s) of a sub-period run
    keep = (steps.sum(axis=0) >= min_side_months) & ((1 - steps).sum(axis=0) >= min_side_months)
    Xs = np.column_stack([X, steps[:, keep]])
    fit_s = sm.OLS(y, Xs).fit()
    se_s = np.sqrt(np.diag(_gap_aware_newey_west(Xs, fit_s.resid, mi, lags[0])))
    step_info = {}
    for i, (s, v) in enumerate(zip(np.array(PROCESSOR_STEPS)[keep], fit_s.params[6:])):
        step_info[f"step_{s}_ppb"] = float(v)
        step_info[f"step_{s}_se_nw"] = float(se_s[6 + i])
    step_info["slope_with_processor_steps_se_nw"] = float(se_s[1])
    step_info["corr_slope_last_step"] = float(np.corrcoef(Xs[:, 1], Xs[:, -1])[0, 1]) if keep.any() else None
    # alternative specification: a change of SLOPE (hinge) at the 2022-07-17
    # bias-correction change instead of level steps. The data cannot separate
    # a level step from a change in growth rate; AIC is reported for each.
    step_info["aic_trend_only"] = float(fit.aic)
    step_info["aic_with_steps"] = float(fit_s.aic)
    hinge_at = pd.Timestamp(PROCESSOR_STEPS[-1])
    th = ((hinge_at - pd.Timestamp("2020-01-01")).days / 365.25)
    after = (d["date"] >= hinge_at).values
    if after.sum() >= min_side_months and (~after).sum() >= min_side_months:
        Xh = np.column_stack([X, np.clip(t - th, 0, None)])
        fit_h = sm.OLS(y, Xh).fit()
        d["fitted_hinge"] = fit_h.fittedvalues
        se_h = np.sqrt(np.diag(_gap_aware_newey_west(Xh, fit_h.resid, mi, lags[0])))
        step_info.update(slope_before_hinge_ppb_per_yr=float(fit_h.params[1]),
                         slope_before_hinge_se_nw=float(se_h[1]),
                         slope_change_at_hinge_ppb_per_yr=float(fit_h.params[6]),
                         slope_change_at_hinge_se_nw=float(se_h[6]),
                         aic_with_hinge=float(fit_h.aic))
    out = dict(slope_ppb_per_yr=float(fit.params[1]), slope_se=se[f"slope_se_nw_lag{lags[0]}"], **se,
               slope_with_processor_steps_ppb_per_yr=float(fit_s.params[1]), **step_info,
               seasonal_amplitude_ppb=float((S.max() - S.min()) / 2),
               annual_harmonic_amplitude_ppb=float(np.hypot(fit.params[2], fit.params[3])),
               semiannual_harmonic_amplitude_ppb=float(np.hypot(fit.params[4], fit.params[5])),
               n_months_used=int(ok.sum()), n_months_dropped=int((~ok).sum()),
               months_dropped=bg.loc[~ok, "month"].tolist(), r2=float(fit.rsquared))
    return out, d


def cluster_monthly_series(labels, monthly_anom: dict, monthly_count: dict, clusters):
    """Count-weighted mean de-trended anomaly over each cluster's pixels, per month."""
    rows = []
    for k in clusters:
        m = labels == k
        for key in sorted(monthly_anom):
            a = monthly_anom[key][m]
            n = monthly_count.get(key, np.ones_like(monthly_anom[key]))[m]
            ok = np.isfinite(a) & np.isfinite(n) & (n > 0)
            val = float(np.sum(a[ok] * n[ok]) / np.sum(n[ok])) if ok.any() else np.nan
            rows.append(dict(cluster=int(k), month=key, anom_ppb=val,
                             n_retrievals=float(np.sum(n[ok])) if ok.any() else 0.0,
                             n_pixels_with_data=int(ok.sum())))
    df = pd.DataFrame(rows)
    df["date"] = pd.to_datetime(df["month"] + "-15")
    return df
