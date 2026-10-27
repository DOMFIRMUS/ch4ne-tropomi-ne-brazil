"""Real-data check (interactive, 1 x 1 deg box in Bahia, 2023): the one-pass
statistics image of ch4ne.gee_export.stats_image against an independent
two-pass computation, for ALL 16 bands that come from the retrievals, and the
effect of countDistinct vs countDistinctNonNull on n_months (cell by cell,
against the monthly stack). No synthetic data. Output: verifica_passo_unico_v2.log"""
import json
import sys
from pathlib import Path as _P

sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
import ee  # noqa: E402

from ch4ne import gee_export as G  # noqa: E402
from ch4ne.config import STATS_CORE_BANDS, Config  # noqa: E402

cfg = Config(base_dir=_P("."))
G.initialize(cfg)
proj = G.grid_projection(cfg)
box = ee.Geometry.Rectangle([-38.5, -12.5, -37.5, -11.5], proj="EPSG:4326", geodesic=False)
kw = dict(geometry=box, crs="EPSG:4326", crsTransform=G.crs_transform(cfg), maxPixels=1e8)
anom = G.anomaly_collection(cfg, G.orbit_collection(cfg), G.background_table(cfg)).filterDate("2023-01-01", "2024-01-01")

# ---- new: one pass (as exported)
new = G.stats_image(cfg, anom)

# ---- reference: separate passes, single-input reducers only
pct = anom.select("anom").reduce(ee.Reducer.percentile([5, 10, 50, 90, 95]).combine(
    ee.Reducer.mean(), sharedInputs=True).combine(ee.Reducer.stdDev(), sharedInputs=True).combine(
    ee.Reducer.count(), sharedInputs=True))
means = anom.select(["xch4", "t", "janjun", "aod"]).reduce(ee.Reducer.mean())
p = {q: pct.select(f"anom_p{q}") for q in (5, 10, 50, 90, 95)}
keys = [f"2023-{m:02d}" for m in range(1, 13)]
stack_n = G.monthly_stack(cfg, anom, keys).select("n_.*").unmask(0).gt(0).reduce(ee.Reducer.sum())
ref = ee.Image.cat([
    pct.select("anom_count"), means.select("xch4_mean"), pct.select("anom_mean"), p[50], pct.select("anom_stdDev"),
    p[10], p[90], p[5], p[95], p[90].subtract(p[50]), p[50].subtract(p[10]), p[95].subtract(p[50]),
    p[50].subtract(p[5]), means.select("t_mean"), means.select("janjun_mean"), means.select("aod_mean"), stack_n,
]).rename(list(STATS_CORE_BANDS)).toFloat()

diff = new.subtract(ref).abs()
# countDistinct (plain) on the same month index, to document the +1 effect cell by cell
plain = anom.select(["anom", "mi"]).reduce(ee.Reducer.count().combine(
    ee.Reducer.countDistinct().setOutputs(["nd"]), sharedInputs=False)).select("nd")
d_plain = plain.subtract(stack_n).rename("d_plain")
d_nonnull = new.select("n_months").subtract(stack_n).rename("d_nonnull")
hist = ee.Image.cat([d_plain, d_nonnull]).updateMask(new.select("n_obs").gt(0))
out = ee.Dictionary({
    "cells_with_data": new.select("n_obs").reduceRegion(ee.Reducer.count(), **kw),
    "maxabs_diff_one_pass_vs_reference": diff.reduceRegion(ee.Reducer.max(), **kw),
    "countDistinct_minus_stack_histogram": hist.select("d_plain").reduceRegion(ee.Reducer.frequencyHistogram(), **kw),
    "countDistinctNonNull_minus_stack_histogram": hist.select("d_nonnull").reduceRegion(
        ee.Reducer.frequencyHistogram(), **kw),
}).getInfo()
text = json.dumps(out, indent=1, sort_keys=True)
print(text)
(_P(__file__).with_suffix(".log")).write_text(text, encoding="utf-8")
