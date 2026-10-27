import sys, json, traceback
from pathlib import Path as _P; sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from pathlib import Path
import ee
from ch4ne.config import Config
from ch4ne import gee_export as G
cfg = Config(base_dir=Path("."))
G.initialize(cfg)
proj = G.grid_projection(cfg)
box = ee.Geometry.Rectangle([-38.0, -9.5, -37.5, -9.0], proj="EPSG:4326", geodesic=False)
kw = dict(reducer=ee.Reducer.mean(), geometry=box, crs="EPSG:4326", crsTransform=G.crs_transform(cfg), maxPixels=1e8)
cov = G.covariates(cfg)
occ = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence")
print("GSW proj:", occ.projection().getInfo())
water = occ.gte(50).unmask(0).rename("w")
print("water proj:", water.projection().getInfo())
w300 = water.reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(crs="EPSG:4326", scale=300)
tests = {
  "albedo": cov.select("albedo_b7"),
  "elev": cov.select(["elev_mean", "elev_sd"]),
  "w300_at_300m": w300,
  "water_frac": cov.select("water_frac"),
  "all_unmasked": cov.unmask(-9999, False),
}
for k, img in tests.items():
    try:
        if k == "w300_at_300m":
            r = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=box, scale=300, maxPixels=1e8).getInfo()
        else:
            r = img.reduceRegion(**kw).getInfo()
        print(k, "OK", r)
    except Exception as e:
        print(k, "ERROR", str(e)[:300])
