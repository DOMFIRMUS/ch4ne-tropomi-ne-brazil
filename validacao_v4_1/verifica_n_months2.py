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
anom = G.anomaly_collection(cfg, G.orbit_collection(cfg), G.background_table(cfg)).filterDate("2023-01-01", "2024-01-01")
s_new = G.stats_image(cfg, anom)
keys = [f"2023-{m:02d}" for m in range(1, 13)]
mon = G.monthly_stack(cfg, anom, keys).select("n_.*").unmask(0)
nm_stack = mon.gt(0).reduce(ee.Reducer.sum()).rename("nm_stack")
img = (s_new.select(["n_months", "n_obs"]).addBands(nm_stack).addBands(mon)
       .addBands(ee.Image.pixelLonLat()))
rows = img.sample(region=box, projection=proj, geometries=False, dropNulls=False).getInfo()["features"]
bad = [r["properties"] for r in rows if r["properties"].get("n_months") is not None
       and r["properties"]["n_months"] != r["properties"]["nm_stack"]]
print("cells sampled:", len(rows), "with n_months:", sum(r["properties"].get("n_months") is not None for r in rows),
      "different:", len(bad))
for p in bad[:4]:
    print(json.dumps(p))
