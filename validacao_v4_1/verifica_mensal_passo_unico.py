"""Real data, small box: one-pass monthly stack vs the previous filtered sub-collection version (2023)."""
import sys, json
from pathlib import Path as _P; sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from pathlib import Path
import ee
from ch4ne.config import Config
from ch4ne import gee_export as G
cfg = Config(base_dir=Path("."))
G.initialize(cfg)
box = ee.Geometry.Rectangle([-38.5, -12.5, -37.5, -11.5], proj="EPSG:4326", geodesic=False)
kw = dict(geometry=box, crs="EPSG:4326", crsTransform=G.crs_transform(cfg), maxPixels=1e8)
anom = G.anomaly_collection(cfg, G.orbit_collection(cfg), G.background_table(cfg))
keys = [f"2023-{m:02d}" for m in range(1, 13)]

def old_stack(anom, keys):
    empty = ee.Image.constant(0).rename("anom").updateMask(0).toFloat()
    means, counts = [], []
    for key in keys:
        sub = anom.filter(ee.Filter.eq("month", key)).select("anom").merge(ee.ImageCollection([empty]))
        tag = key.replace("-", "_")
        r = sub.reduce(ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True))
        means.append(r.select("anom_mean").rename(f"a_{tag}"))
        counts.append(r.select("anom_count").rename(f"n_{tag}"))
    return ee.Image.cat(means + counts).toFloat()

new = G.monthly_stack(cfg, anom, keys)
old = old_stack(anom, keys)
print("bands equal:", new.bandNames().getInfo() == old.bandNames().getInfo())
d = new.unmask(-1e6).subtract(old.unmask(-1e6)).abs()
out = ee.Dictionary({
    "maxabs_diff": d.reduceRegion(ee.Reducer.max(), **kw),
    "valid_new": new.select("a_.*").reduceRegion(ee.Reducer.count(), **kw),
    "valid_old": old.select("a_.*").reduceRegion(ee.Reducer.count(), **kw),
    "total_count_new": new.select("n_.*").reduceRegion(ee.Reducer.sum(), **kw),
}).getInfo()
print("max |new-old| over all 24 bands:", max(v for v in out["maxabs_diff"].values()))
print("valid cells per month new == old:", out["valid_new"] == out["valid_old"])
print(json.dumps({k: out["valid_new"][k] for k in sorted(out["valid_new"])}))
print("retrievals per month (sum over box):", {k[2:]: round(v) for k, v in sorted(out["total_count_new"].items())})
