"""
Provenance of the Earth Engine exports.

The analysis only runs on rasters produced by gee_export.py from the real
COPERNICUS/S5P/OFFL/L3_CH4 collection: download_all() writes PROVENANCE.json
(asset ids and update times, the run identity stamped on the assets,
processor versions, band names and the SHA-256 of every file, band-name
sidecars included). analysis.run() refuses to start if that record is
missing, if any file changed after download, or if the record was produced
with settings that differ from the configuration being analysed.
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

RECORD = "PROVENANCE.json"


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_record(gee_dir, files: dict, assets: dict, cfg, extra=None):
    import ee
    rec = dict(
        created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        source="Google Earth Engine",
        gee_project=cfg.gee_project,
        collection=cfg.collection,
        band=cfg.band,
        start_date=cfg.start_date,
        end_date=cfg.end_date,
        run_label=cfg.run_label,
        earthengine_api=ee.__version__,
        assets=assets,
        files={name: dict(path=Path(p).name, sha256=sha256(p)) for name, p in files.items()},
        config=cfg.as_dict(),
    )
    rec.update(extra or {})
    Path(gee_dir, RECORD).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    return rec


def verify(gee_dir, cfg=None) -> dict:
    """Raise unless every exported file is present, unchanged since download
    and (when cfg is given) produced with the same run identity and settings."""
    gee_dir = Path(gee_dir)
    p = gee_dir / RECORD
    if not p.exists():
        raise RuntimeError(
            f"{p} not found. The analysis only runs on rasters downloaded from Earth Engine "
            "by `python scripts/run_pipeline.py download` (which writes this provenance record).")
    rec = json.loads(p.read_text(encoding="utf-8"))
    if rec.get("source") != "Google Earth Engine" or not rec.get("assets"):
        raise RuntimeError(f"{p} does not describe an Earth Engine export.")
    for name, meta in rec["files"].items():
        f = gee_dir / meta["path"]
        if not f.exists():
            raise RuntimeError(f"Exported file missing: {f}")
        if sha256(f) != meta["sha256"]:
            raise RuntimeError(f"{f} differs from the file downloaded from Earth Engine (checksum mismatch).")
    if cfg is not None:
        mine = (cfg.start_date, cfg.end_date, cfg.run_label)
        theirs = (rec.get("start_date"), rec.get("end_date"), rec.get("run_label", ""))
        if mine != theirs:
            raise RuntimeError(f"{p} is for period/label {theirs}, but the analysis is configured for {mine}.")
        params = rec.get("run_properties", {}).get("ch4ne_params")
        if params is not None and json.loads(params) != json.loads(json.dumps(cfg.gee_params(), sort_keys=True)):
            raise RuntimeError(f"{p} was exported with different Earth Engine settings than ch4ne/config.py now "
                               "has; re-export or restore the settings.")
    return rec
