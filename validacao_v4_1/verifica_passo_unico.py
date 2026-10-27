"""Real-data check (interactive, small box): new single-reduceResolution orbit prep and
single-pass stats reducer vs the previous two-step versions. No synthetic data."""
import sys, json
from pathlib import Path as _P; sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from pathlib import Path
import ee
from ch4ne.config import Config
from ch4ne import gee_export as G

cfg = Config(base_dir=Path("."))
G.initialize(cfg)
proj = G.grid_projection(cfg)
box = ee.Geometry.Rectangle([-38.5, -12.5, -37.5, -11.5], proj="EPSG:4326", geodesic=False)
kw = dict(geometry=box, crs="EPSG:4326", crsTransform=G.crs_transform(cfg), maxPixels=1e8)

def old_prep(img):
    ok = (img.select(cfg.uncertainty_band).lte(cfg.max_uncertainty_ppb)
          .And(img.select(cfg.aod_band).lte(cfg.max_aod_swir)))
    ch4 = img.select(cfg.band).updateMask(ok)
    aod = img.select(cfg.aod_band).updateMask(ch4.mask())
    agg = (ch4.addBands(aod).rename(["xch4", "aod"]).reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(proj))
    frac = (ch4.mask().gt(0).rename("f").reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(proj))
    return agg.updateMask(frac.gte(cfg.min_cell_valid_fraction))

# ---- test 1: orbit prep on 6 real orbits crossing the box
src = G.source_collection(cfg, "2023-08-01", "2023-08-08").filterBounds(box)
new = G.orbit_collection(cfg).filterDate("2023-08-01", "2023-08-08").filterBounds(box)
ids = src.aggregate_array("system:index").slice(0, 6)
res = []
for i in ids.getInfo():
    o = old_prep(src.filter(ee.Filter.eq("system:index", i)).first())
    n = ee.Image(new.filter(ee.Filter.eq("system:index", i)).first())
    d = n.select(["xch4", "aod"]).subtract(o).abs().rename(["dx", "da"])
    st = ee.Dictionary({
        "n_old": o.select("xch4").reduceRegion(ee.Reducer.count(), **kw).get("xch4"),
        "n_new": n.select("xch4").reduceRegion(ee.Reducer.count(), **kw).get("xch4"),
        "maxdiff": d.reduceRegion(ee.Reducer.max(), **kw)}).getInfo()
    res.append((i, st))
    print("orbit", i, st, flush=True)

# ---- test 2: stats single pass vs separate passes (2023, same box)
c23 = Config(base_dir=Path("."))
anom = G.anomaly_collection(cfg, G.orbit_collection(cfg), G.background_table(cfg)).filterDate("2023-01-01", "2024-01-01")
s_new = G.stats_image(cfg, anom)
print("new bands:", s_new.bandNames().getInfo())
red = (ee.Reducer.mean().combine(ee.Reducer.stdDev(), sharedInputs=True).combine(ee.Reducer.count(), sharedInputs=True)
       .combine(ee.Reducer.percentile([5, 10, 50, 90, 95], maxRaw=cfg.percentile_max_raw), sharedInputs=True))
s_old = anom.select("anom").reduce(red, cfg.parallel_scale)
m_old = anom.select(["xch4", "t", "janjun", "aod"]).reduce(ee.Reducer.mean(), cfg.parallel_scale)
old = ee.Image.cat([s_old.select("anom_count"), m_old.select("xch4_mean"), s_old.select("anom_mean"),
                    s_old.select("anom_p50"), s_old.select("anom_stdDev"), s_old.select("anom_p10"), s_old.select("anom_p90"),
                    s_old.select("anom_p5"), s_old.select("anom_p95")]).rename(
    ["n_obs", "xch4_mean", "anom_mean", "anom_median", "anom_sd", "anom_p10", "anom_p90", "anom_p05", "anom_p95"])
mon = G.monthly_stack(cfg, anom, [f"2023-{m:02d}" for m in range(1, 13)])
nm_stack = mon.select("n_.*").gt(0).reduce(ee.Reducer.sum()).rename("nm_stack")
cmp_b = old.bandNames()
diff = s_new.select(cmp_b).subtract(old).abs()
out = ee.Dictionary({
    "maxabs_diff_by_band": diff.reduceRegion(ee.Reducer.max(), **kw),
    "cells_new": s_new.select("n_obs").reduceRegion(ee.Reducer.count(), **kw),
    "cells_old": old.select("n_obs").reduceRegion(ee.Reducer.count(), **kw),
    "n_months_minus_stack_maxabs": s_new.select("n_months").subtract(nm_stack).abs().reduceRegion(ee.Reducer.max(), **kw),
    "n_months_range": s_new.select("n_months").reduceRegion(ee.Reducer.minMax(), **kw),
    "n_obs_range": s_new.select("n_obs").reduceRegion(ee.Reducer.minMax(), **kw),
    "anom_mean_range": s_new.select("anom_mean").reduceRegion(ee.Reducer.minMax(), **kw),
}).getInfo()
print(json.dumps(out, indent=1))
