"""
Publication figures (MDPI Remote Sensing requirements).

* No titles or figure numbers inside the images (captions live in the paper).
* PNG at cfg.dpi (600) plus vector PDF; sized for 17.5 cm (full) text width.
* Every map: graticule labelled in degrees S/W, geodesic scale bar, north
  arrow, equal-area-looking aspect (1/cos(lat)), rasters drawn pixel-exact
  from the GeoTIFF transform (v3 drew Earth Engine thumbnails with bilinear
  smoothing on a lat/lon extent).
* Perceptually uniform / colour-blind-safe palettes (no 'jet').
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MultipleLocator
from pyproj import Geod

GEOD = Geod(ellps="WGS84")
# Okabe-Ito palette
OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000", "grey": "#999999"}
CLASS_COLORS = {"Oil pipeline": OI["vermillion"], "Gathering gas pipeline": OI["orange"],
                "Transport gas pipeline": OI["blue"], "Distribution gas pipeline": OI["green"]}


def setup_style(cfg):
    mpl.rcParams.update({
        "font.size": cfg.font_size, "axes.titlesize": cfg.font_size + 1, "axes.labelsize": cfg.font_size,
        "xtick.labelsize": cfg.font_size - 1, "ytick.labelsize": cfg.font_size - 1,
        "legend.fontsize": cfg.font_size - 1, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "savefig.dpi": cfg.dpi, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02, "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    })


def save(fig, out_dir, stem, cfg):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.png", dpi=cfg.dpi)
    fig.savefig(out_dir / f"{stem}.pdf")
    plt.close(fig)
    return out_dir / f"{stem}.png"


# ------------------------------------------------------------------ map helpers
def _fmt_lon(x, _):
    return f"{abs(x):g}°{'W' if x < 0 else ('E' if x > 0 else '')}"


def _fmt_lat(y, _):
    return f"{abs(y):g}°{'S' if y < 0 else ('N' if y > 0 else '')}"


def map_frame(ax, bounds, step=4):
    """bounds = (west, south, east, north)."""
    w, s, e, n = bounds
    ax.set_xlim(w, e)
    ax.set_ylim(s, n)
    ax.set_aspect(1.0 / np.cos(np.deg2rad((s + n) / 2)))
    ax.xaxis.set_major_locator(MultipleLocator(step))
    ax.yaxis.set_major_locator(MultipleLocator(step))
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_lon))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_lat))
    ax.grid(True, lw=0.3, ls=":", color="0.5", zorder=0)
    ax.tick_params(direction="out", length=2)


def scale_bar(ax, length_km=200, loc=(0.60, 0.05), segments=2, lw=3):
    """Geodesic scale bar drawn along the parallel at the bar's position."""
    x0f, y0f = loc
    (w, e), (s, n) = ax.get_xlim(), ax.get_ylim()
    lon0, lat0 = w + x0f * (e - w), s + y0f * (n - s)
    step = length_km / segments
    xs = [lon0]
    for i in range(segments):
        lon1, _, _ = GEOD.fwd(xs[-1], lat0, 90.0, step * 1000.0)
        xs.append(lon1)
    for i in range(segments):
        ax.plot(xs[i:i + 2], [lat0, lat0], color="k" if i % 2 == 0 else "w", lw=lw,
                solid_capstyle="butt", zorder=30)
    ax.plot([xs[0], xs[-1]], [lat0, lat0], color="k", lw=lw + 1.2, solid_capstyle="butt", zorder=29)
    for x, lab, ha in ((xs[0], "0", "center"), (xs[-1], f"{length_km:g} km", "center")):
        ax.text(x, lat0 + 0.012 * (n - s), lab, ha=ha, va="bottom", fontsize=mpl.rcParams["font.size"] - 1.5,
                zorder=30)


def north_arrow(ax, loc=(0.92, 0.88), size=0.07):
    x, y = loc
    ax.annotate("N", xy=(x, y + size), xytext=(x, y), xycoords="axes fraction", ha="center", va="center",
                fontsize=mpl.rcParams["font.size"] + 1, fontweight="bold", zorder=30,
                arrowprops=dict(facecolor="k", edgecolor="k", width=1.5, headwidth=6, headlength=6))


def raster(ax, arr, grid, **kw):
    kw.setdefault("interpolation", "nearest")
    im = ax.imshow(np.ma.masked_invalid(arr), extent=grid.extent, origin="upper", aspect="auto", zorder=2, **kw)
    return im


def ne_bounds(ne, pad=0.3):
    """Frame of the NE mainland (islands such as Fernando de Noronha, 32.4 W,
    would otherwise widen every map by ~15 % of empty ocean)."""
    import shapely
    main = ne.geometry.map(lambda g: max(shapely.get_parts(g), key=lambda p: p.area))
    w, s, e, n = main.total_bounds
    return (w - pad, s - pad, e + pad, n + pad)


def _base(ax, ne, lines=None, fill="#f2f2f2", lines_kw=None):
    ne.plot(ax=ax, facecolor=fill, edgecolor="none", zorder=1)
    ne.boundary.plot(ax=ax, color="0.25", lw=0.4, zorder=5)
    if lines is not None:
        kw = dict(color="k", lw=0.5, zorder=6)
        kw.update(lines_kw or {})
        lines.plot(ax=ax, **kw)


def _cbar(fig, ax, im, label, **kw):
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.045, pad=0.07, **kw)
    cb.set_label(label)
    cb.outline.set_linewidth(0.5)
    return cb


def _panel(ax, letter):
    """Panel letter above the top-left corner (never over data, legends or ticks)."""
    ax.set_title(f"({letter})", loc="left", fontweight="bold", fontsize=mpl.rcParams["font.size"], pad=3)


# ------------------------------------------------------------------- figures
def fig_study_area(ne, states, infra, lines, cfg, out_dir, inset=True):
    """inset=False drops the locator map of Brazil (Elsevier maps may not show more than the study area)."""
    fig, ax = plt.subplots(figsize=(cfg.fig_width_full_in * 0.62, cfg.fig_width_full_in * 0.78))
    b = ne_bounds(ne)
    states.plot(ax=ax, facecolor="#fafafa", edgecolor="0.7", lw=0.3, zorder=0)
    _base(ax, ne, fill="#ececec")
    other = infra[infra.location != "onshore NE"]
    if len(other):
        other.plot(ax=ax, color="0.55", lw=0.6, linestyle="--", zorder=6)
    handles = []
    for cls, col in CLASS_COLORS.items():
        sub = lines[lines.asset_class == cls]
        if len(sub):
            sub.plot(ax=ax, color=col, lw=0.9, zorder=7)
            n_seg = int(((infra.asset_class == cls) & (infra.location == "onshore NE")).sum())
            handles.append(Line2D([0], [0], color=col, lw=1.5, label=f"{cls} (n = {n_seg})"))
    handles.append(Line2D([0], [0], color="0.55", lw=1, ls="--", label="Offshore / outside NE"))
    for _, s in ne.iterrows():
        p = s.geometry.representative_point()
        ax.text(p.x, p.y, s.UF, ha="center", va="center", fontsize=cfg.font_size - 1, color="0.3", zorder=8)
    map_frame(ax, b)
    scale_bar(ax, 200, loc=(0.05, 0.04))
    north_arrow(ax)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, frameon=False,
              handlelength=1.5, columnspacing=1.0)
    if inset:
        # locator inset over the Atlantic (lower right)
        ins = ax.inset_axes([0.73, 0.02, 0.26, 0.26])
        states.plot(ax=ins, facecolor="0.85", edgecolor="0.6", lw=0.2)
        ne.plot(ax=ins, facecolor="0.35", edgecolor="none")
        ins.set_xticks([]); ins.set_yticks([])
        ins.set_aspect(1.0)
        for sp in ins.spines.values():
            sp.set_linewidth(0.5)
    return save(fig, out_dir, "Fig01_study_area_infrastructure", cfg)


def fig_composite_and_density(f, grid, ne, lines, valid_region, cfg, out_dir):
    fig, axs = plt.subplots(1, 2, figsize=(cfg.fig_width_full_in, cfg.fig_width_full_in * 0.66))
    b = ne_bounds(ne)
    x = np.where(valid_region, f["xch4_mean"], np.nan)
    lo, hi = np.nanpercentile(x, [1, 99])
    for ax, arr, cmap, lab, norm, letter in (
            (axs[0], x, "cividis", f"XCH$_4$ mean, {cfg.period_label} (ppb)", mpl.colors.Normalize(lo, hi), "a"),
            (axs[1], np.where(valid_region, f["n_obs"], np.nan), "magma",
             "Valid retrievals per cell", mpl.colors.LogNorm(max(1, np.nanpercentile(f["n_obs"], 1)),
                                                             np.nanpercentile(f["n_obs"], 99.5)), "b")):
        _base(ax, ne, lines, lines_kw=dict(color="w", lw=0.5))
        im = raster(ax, arr, grid, cmap=cmap, norm=norm)
        map_frame(ax, b)
        scale_bar(ax, 200, loc=(0.05, 0.04))
        north_arrow(ax)
        _panel(ax, letter)
        _cbar(fig, ax, im, lab, extend="both" if letter == "a" else "neither")
    return save(fig, out_dir, "Fig02_composite_and_data_density", cfg)


def fig_background(bg_fit_df, trend, bg_all, cfg, out_dir):
    fig, ax = plt.subplots(figsize=(cfg.fig_width_full_in, 3.3))
    ax.plot(bg_all["date"], bg_all["bg_ppb"], "o", ms=2.5, color="0.6", label="Month dropped (coverage)")
    ax.plot(bg_fit_df["date"], bg_fit_df["bg_ppb"], "o", ms=2.5, color=OI["blue"], label="Monthly NE background")
    ax.plot(bg_fit_df["date"], bg_fit_df["fitted"], "-", lw=1, color=OI["vermillion"],
            label=f"Linear trend + seasonal cycle ({trend['slope_ppb_per_yr']:.1f} ppb yr$^{{-1}}$)")
    if "fitted_hinge" in bg_fit_df and "slope_before_hinge_ppb_per_yr" in trend:
        a = trend["slope_before_hinge_ppb_per_yr"]
        b = a + trend["slope_change_at_hinge_ppb_per_yr"]
        ax.plot(bg_fit_df["date"], bg_fit_df["fitted_hinge"], "--", lw=1, color="k",
                label=f"Slope change at 2022-07-17 ({a:.1f} → {b:.1f} ppb yr$^{{-1}}$)")
    ax.axvline(pd.Timestamp("2022-07-17"), color="0.4", lw=0.8, ls=":")
    ax.text(pd.Timestamp("2022-07-17"), 0.03, " TROPOMI bias-correction change", transform=ax.get_xaxis_transform(),
            fontsize=cfg.font_size - 2, color="0.3", va="bottom")
    ax.set_ylabel("XCH$_4$ (ppb)")
    ax.grid(True, lw=0.3, ls=":")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=cfg.font_size - 1)
    fig.tight_layout()
    return save(fig, out_dir, "Fig03_regional_background_trend", cfg)


def fig_anomaly_maps(f, grid, ne, lines, polys, screened, cfg, out_dir, z_thr):
    fig, axs = plt.subplots(1, 2, figsize=(cfg.fig_width_full_in, cfg.fig_width_full_in * 0.66))
    b = ne_bounds(ne)
    a = np.where(screened, f["anom_mean"], np.nan)
    lim = np.nanpercentile(np.abs(a - np.nanmedian(a)), 99)
    med = np.nanmedian(a)
    excluded = np.where(np.isfinite(f["anom_mean"]) & ~screened & f["in_ne"], 1.0, np.nan)
    for ax, arr, cmap, norm, lab, letter in (
            (axs[0], a, "RdBu_r", mpl.colors.TwoSlopeNorm(med, med - lim, med + lim),
             "De-trended anomaly (ppb)", "a"),
            (axs[1], np.where(screened, f["z"], np.nan), "RdBu_r", mpl.colors.TwoSlopeNorm(0, -4, 4),
             "Robust z of local enhancement", "b")):
        _base(ax, ne, None)
        raster(ax, excluded, grid, cmap=mpl.colors.ListedColormap(["0.75"]), vmin=0, vmax=1)
        im = raster(ax, arr, grid, cmap=cmap, norm=norm)
        lines.plot(ax=ax, color="k", lw=0.45, zorder=6)
        if letter == "b" and polys is not None and len(polys):
            polys.boundary.plot(ax=ax, color=OI["green"], lw=0.8, zorder=7)
        map_frame(ax, b)
        scale_bar(ax, 200, loc=(0.05, 0.04))
        north_arrow(ax)
        _panel(ax, letter)
        _cbar(fig, ax, im, lab, extend="both")
    fig.legend(handles=[Patch(fc="0.75", label="Screened out"),
                        Line2D([0], [0], color=OI["green"], lw=1, label=f"Cluster (z > {z_thr:g})"),
                        Line2D([0], [0], color="k", lw=0.8, label="Onshore pipeline")],
               loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=False,
               fontsize=cfg.font_size - 1)
    return save(fig, out_dir, "Fig04_anomaly_and_clusters", cfg)


def fig_histograms(legacy_x, legacy_mu, legacy_sd, enh, enh_med, enh_mad, z_thr, cfg, out_dir):
    fig, axs = plt.subplots(1, 2, figsize=(cfg.fig_width_full_in, 3.4))
    ax = axs[0]
    x = legacy_x[np.isfinite(legacy_x)]
    ax.hist(x, bins=80, color=OI["sky"], edgecolor="none")
    ax.axvline(legacy_mu, color="k", ls="--", lw=0.8, label=f"Mean = {legacy_mu:.2f} ppb")
    ax.axvline(legacy_mu + 2 * legacy_sd, color=OI["vermillion"], lw=1,
               label=f"Mean + 2σ = {legacy_mu + 2 * legacy_sd:.2f} ppb")
    ax.set_yscale("log")
    ax.set_xlabel("XCH$_4$ composite (ppb)")
    ax.set_ylabel(f"Number of cells (n = {x.size:,})")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=cfg.font_size - 1)
    _panel(ax, "a")
    ax = axs[1]
    e = enh[np.isfinite(enh)]
    ax.hist(e, bins=80, color=OI["green"], edgecolor="none")
    ax.axvline(enh_med, color="k", ls="--", lw=0.8, label=f"Median = {enh_med:.2f} ppb")
    ax.axvline(enh_med + z_thr * enh_mad, color=OI["vermillion"], lw=1,
               label=f"Median + {z_thr:g}·1.4826·MAD = {enh_med + z_thr * enh_mad:.2f} ppb")
    ax.set_yscale("log")
    ax.set_xlabel("Enhancement over local background (ppb)")
    ax.set_ylabel(f"Number of cells (n = {e.size:,})")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=cfg.font_size - 1)
    _panel(ax, "b")
    fig.subplots_adjust(wspace=0.35)
    return save(fig, out_dir, "Fig05_distributions", cfg)


def fig_states(df_pix, cfg, out_dir, var="xch4_mean", ylabel=None, stem="Fig06_states"):
    ylabel = ylabel or f"XCH$_4$ mean, {cfg.period_label} (ppb)"
    order = df_pix.groupby("UF")[var].median().sort_values(ascending=False).index.tolist()
    data = [df_pix.loc[df_pix.UF == u, var].dropna().values for u in order]
    fig, ax = plt.subplots(figsize=(cfg.fig_width_full_in, 2.6))
    parts = ax.violinplot(data, showextrema=False, widths=0.85)
    for pc in parts["bodies"]:
        pc.set_facecolor(OI["sky"]); pc.set_edgecolor("0.3"); pc.set_linewidth(0.4); pc.set_alpha(0.8)
    ax.boxplot(data, widths=0.18, showfliers=False, medianprops=dict(color="k", lw=1),
               boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.axhline(np.nanmean(df_pix[var]), color="k", ls="--", lw=0.7, label="Northeast mean")
    ax.set_xticks(range(1, len(order) + 1), [f"{u}\n{len(d):,}" for u, d in zip(order, data)],
                  fontsize=cfg.font_size - 2)
    ax.set_xlabel("State (number of cells)")
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y", lw=0.3, ls=":")
    ax.legend(frameon=False, loc="lower left")
    return save(fig, out_dir, stem, cfg)


def fig_rings(ring_tab, order, cfg, out_dir, panels):
    """panels: list of (variable, label). Median with IQR per ring."""
    n = len(panels)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(cfg.fig_width_full_in, 2.0 * nrow), squeeze=False)
    xs = np.arange(len(order))
    for i, (var, lab) in enumerate(panels):
        ax = axs.flat[i]
        t = ring_tab[ring_tab.variable == var].set_index("ring").reindex(order)
        ax.errorbar(xs, t["median"], yerr=[t["median"] - t["q25"], t["q75"] - t["median"]],
                    fmt="o", color=OI["blue"], ms=3, lw=0.8, capsize=2, label="Median (IQR)")
        ax.plot(xs, t["mean"], "s", color=OI["vermillion"], ms=2.5, label="Mean")
        ax.set_xticks(xs, [o.replace(" km", "") for o in order], rotation=0)
        ax.set_xlabel("Distance to pipeline (km)")
        ax.set_ylabel(lab)
        ax.grid(True, lw=0.3, ls=":")
        _panel(ax, chr(97 + i))
    for j in range(n, nrow * ncol):
        axs.flat[j].axis("off")
    h, lab = axs.flat[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, fontsize=cfg.font_size - 1)
    return save(fig, out_dir, "Fig07_proximity_rings", cfg)


def fig_cluster_timeseries(ts, clusters, cfg, out_dir, n_show=6):
    ks = list(clusters)[:n_show]
    if not ks:
        return None
    ncol = 2
    nrow = int(np.ceil(len(ks) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(cfg.fig_width_full_in, 1.6 * nrow), sharex=True, squeeze=False)
    for ax, k in zip(axs.flat, ks):
        d = ts[ts.cluster == k]
        sz = 4 + 20 * (d["n_retrievals"] / max(1.0, d["n_retrievals"].max()))
        ax.scatter(d["date"], d["anom_ppb"], s=sz, color=OI["blue"], lw=0, zorder=3)
        ax.plot(d["date"], d["anom_ppb"].rolling(6, center=True, min_periods=3).mean(), color=OI["vermillion"], lw=1)
        ax.axhline(0, color="k", lw=0.5, ls="--")
        ax.set_ylabel("Anomaly (ppb)")
        ax.set_title(f"Cluster {k}", loc="left", fontweight="bold", fontsize=mpl.rcParams["font.size"], pad=3)
        ax.grid(True, lw=0.3, ls=":")
    for j in range(len(ks), nrow * ncol):
        axs.flat[j].axis("off")
    fig.tight_layout()
    return save(fig, out_dir, "Fig08_cluster_timeseries", cfg)


def fig_sensitivity(sens, cfg, out_dir):
    """Heat maps of the number of clusters (all, and with q below the FDR
    threshold) across data-count rules and z thresholds."""
    if sens is None or sens.empty:
        return None
    sens = sens.assign(rule=[f"obs>={int(a)}, months>={int(b)}" for a, b in zip(sens.min_obs, sens.min_months)])
    order = list(dict.fromkeys(sens["rule"]))
    masks = sens["coastal_screen"].unique()
    fig, axs = plt.subplots(len(masks), 2, figsize=(cfg.fig_width_full_in, 2.4 * len(masks)), squeeze=False)
    for r, cm in enumerate(masks):
        s = sens[sens.coastal_screen == cm]
        for c, (col, lab) in enumerate((("n_clusters", "Clusters"),
                                        ("n_clusters_q", f"Clusters with q <= {cfg.q_threshold:g}"))):
            p = s.pivot(index="rule", columns="z_threshold", values=col).reindex(order)
            ax = axs[r, c]
            vmax = max(1.0, float(np.nanmax(p.values))) if np.isfinite(np.nanmax(p.values)) else 1.0
            im = ax.imshow(p.values, cmap="viridis", aspect="auto", origin="lower", vmin=0, vmax=vmax)
            for (i, j), v in np.ndenumerate(p.values):
                r_, g_, b_, _ = im.cmap(im.norm(v))
                ink = "k" if 0.299 * r_ + 0.587 * g_ + 0.114 * b_ > 0.5 else "w"   # readable on any cell
                ax.text(j, i, f"{int(v)}", ha="center", va="center", color=ink, fontsize=cfg.font_size - 1)
            ax.set_xticks(range(p.shape[1]), [f"{x:g}" for x in p.columns])
            ax.set_yticks(range(p.shape[0]), p.index, fontsize=cfg.font_size - 2)
            ax.set_xlabel("Robust z threshold")
            ax.set_title(f"{lab} (coastal screen {'on' if cm else 'off'})", fontsize=cfg.font_size)
    fig.tight_layout()
    return save(fig, out_dir, "FigS1_sensitivity", cfg)
