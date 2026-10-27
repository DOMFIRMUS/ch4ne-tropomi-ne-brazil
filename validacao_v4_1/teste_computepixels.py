import sys, time, io
from pathlib import Path as _P; sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from pathlib import Path
import numpy as np
import ee
from ch4ne.config import Config
from ch4ne import gee_export as G
cfg = Config(base_dir=Path("."))
G.initialize(cfg)
img = G.covariates(cfg).unmask(-9999, False)
for n, (c0, r0) in ((32, (160, 128)), (64, (160, 128))):   # block over inland Bahia/Pernambuco
    x0 = cfg.grid_west + c0 * cfg.grid_res_deg
    y0 = cfg.grid_north - r0 * cfg.grid_res_deg
    t = time.time()
    arr = ee.data.computePixels({"expression": img, "fileFormat": "NUMPY_NDARRAY",
        "grid": {"dimensions": {"width": n, "height": n},
                 "affineTransform": {"scaleX": cfg.grid_res_deg, "shearX": 0, "translateX": x0,
                                     "shearY": 0, "scaleY": -cfg.grid_res_deg, "translateY": y0},
                 "crsCode": "EPSG:4326"}})
    dt = time.time() - t
    print(n, "px block at", (x0, y0), f"{dt:.1f}s", arr.dtype.names,
          {k: (float(np.nanmin(arr[k])), float(np.nanmax(arr[k]))) for k in arr.dtype.names})
