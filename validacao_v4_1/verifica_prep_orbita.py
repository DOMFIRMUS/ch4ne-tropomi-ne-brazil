import sys, json
from pathlib import Path as _P; sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from pathlib import Path
import ee
from ch4ne.config import Config
from ch4ne import gee_export as G
cfg = Config(base_dir=Path("."))
G.initialize(cfg)
proj = G.grid_projection(cfg)
box = ee.Geometry.Rectangle([-41.0, -11.0, -38.0, -8.0], proj="EPSG:4326", geodesic=False)
kw = dict(geometry=box, crs="EPSG:4326", crsTransform=G.crs_transform(cfg), maxPixels=1e8)

def old_prep(img):
    ok = (img.select(cfg.uncertainty_band).lte(cfg.max_uncertainty_ppb)
          .And(img.select(cfg.aod_band).lte(cfg.max_aod_swir)))
    ch4 = img.select(cfg.band).updateMask(ok)
    aod = img.select(cfg.aod_band).updateMask(ch4.mask())
    agg = (ch4.addBands(aod).rename(["xch4", "aod"]).reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(proj))
    frac = (ch4.mask().gt(0).rename("f").reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(proj))
    return agg.updateMask(frac.gte(cfg.min_cell_valid_fraction))

src = G.source_collection(cfg, "2023-09-01", "2023-09-06").filterBounds(box)
cnt = src.map(lambda i: i.set("n", i.select(cfg.band).reduceRegion(ee.Reducer.count(), box, scale=10000).get(cfg.band)))
ids = cnt.filter(ee.Filter.gt("n", 100)).aggregate_array("system:index").getInfo()[:4]
new = G.orbit_collection(cfg)
for i in ids:
    o = old_prep(src.filter(ee.Filter.eq("system:index", i)).first())
    n = ee.Image(new.filter(ee.Filter.eq("system:index", i)).first())
    d = n.select(["xch4", "aod"]).subtract(o).abs().rename(["dx", "da"])
    st = ee.Dictionary({
        "cells_old": o.select("xch4").reduceRegion(ee.Reducer.count(), **kw).get("xch4"),
        "cells_new": n.select("xch4").reduceRegion(ee.Reducer.count(), **kw).get("xch4"),
        "only_one_side": o.select("xch4").mask().neq(n.select("xch4").mask()).reduceRegion(ee.Reducer.sum(), **kw).get("xch4"),
        "maxabs_diff": d.reduceRegion(ee.Reducer.max(), **kw)}).getInfo()
    print("orbit", i, json.dumps(st), flush=True)
