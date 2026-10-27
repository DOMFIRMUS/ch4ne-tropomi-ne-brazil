"""
Earth Engine stage: the only part of the pipeline that runs on Google servers.

It produces a few small rasters/tables on a fixed 0.0625 deg grid; all
statistics are then computed locally on the full pixel census.

Per orbit (one L3 image = one orbit granule; GEE keeps only qa_value > 0.5):
  1. keep retrievals with uncertainty <= max_uncertainty_ppb and SWIR AOD <=
     max_aod_swir; skip orbits lacking any required band;
  2. AGGREGATE the 0.01 deg L3 bins to the analysis grid (reduceResolution
     mean) and keep a cell only if >= min_cell_valid_fraction of its bins are
     valid (v3 used reproject() alone: nearest-neighbour sampling from an
     image-pyramid level, not an area mean);
Per month (phase A, exported once as a table asset):
  3. regional background B_m = median over NE land of the monthly mean,
     read from the mean pyramid at the grid scale (unweighted, pixel-centre
     rule, exact median);
Per orbit again (phase B, reading the phase-A asset):
  4. anomaly = XCH4 - B_(month of the orbit). This removes the growth and the
     seasonal cycle, so cells observed mostly early/late in the period, or
     mostly in one season, are no longer biased;
Per cell:
  5. mean/median/temporal SD/percentiles/tail widths of the anomaly, number
     of retrievals and of months with data, sampling diagnostics - all in ONE
     pass over the orbits (asset 'stats');
  6. static surface covariates, computed directly on the grid in blocks with
     ee.data.computePixels (no orbit is read; as a batch export the 30 m JRC
     water layer was requested for the whole region at once and failed the
     2^31-pixel limit even with 32-pixel tiles, measured 2026-09-30);
  7. monthly stack (mean anomaly and count per month), one asset per year
     ('monthly_YYYY'), each in one pass over that year's orbits.

Outputs are Earth Engine ASSETS written by batch tasks (no time-out). Every
asset carries the run identity (dates, label, parameters) as properties;
download_all() refuses assets whose identity differs from the configuration
and fetches the stored pixels with ee.data.getPixels on the asset's own grid.
Exports are resumable: an asset that already exists with the same identity,
or a task of the same name still running, is reused instead of recomputed.
"""
import io
import json
import time
import zipfile
from datetime import date, timedelta
from pathlib import Path

import ee

from . import __version__
from .config import COVARIATE_BANDS, NODATA, STATS_BANDS, STATS_CORE_BANDS, Config
from .provenance import write_record

TABLE_JOB = "monthly_background"
ACTIVE_STATES = ("PENDING", "RUNNING")   # a CANCELLING task will not produce the asset


# ------------------------------------------------------------------ basics
def initialize(cfg: Config):
    try:
        ee.Initialize(project=cfg.gee_project)
    except Exception as e:  # not authenticated, or project not enabled
        raise SystemExit(
            f"Earth Engine could not be initialised for project '{cfg.gee_project}': {e}\n"
            "Authenticate once on this machine with:\n"
            "  python -c \"import ee; ee.Authenticate(auth_mode='localhost', force=True)\"\n"
            "and make sure the project is registered at https://code.earthengine.google.com/register")


def crs_transform(cfg: Config):
    return [cfg.grid_res_deg, 0, cfg.grid_west, 0, -cfg.grid_res_deg, cfg.grid_north]


def grid_shape(cfg: Config):
    w = int(round((cfg.grid_east - cfg.grid_west) / cfg.grid_res_deg))
    h = int(round((cfg.grid_north - cfg.grid_south) / cfg.grid_res_deg))
    return w, h


def grid_projection(cfg: Config):
    return ee.Projection("EPSG:4326", crs_transform(cfg))


def export_region(cfg: Config):
    return ee.Geometry.Rectangle([cfg.grid_west, cfg.grid_south, cfg.grid_east, cfg.grid_north],
                                 proj="EPSG:4326", geodesic=False)


def end_exclusive(cfg: Config) -> str:
    return (date.fromisoformat(cfg.end_date) + timedelta(days=1)).isoformat()


def month_keys(cfg: Config):
    y, m = map(int, cfg.start_date[:7].split("-"))
    ye, me = map(int, cfg.end_date[:7].split("-"))
    keys = []
    while (y, m) <= (ye, me):
        keys.append(f"{y:04d}-{m:02d}")
        m += 1
        if m == 13:
            y, m = y + 1, 1
    return keys


def ne_geometry(ne_gdf, tolerance=0.01):
    """NE land as a planar ee.Geometry (pixel-centre inclusion like the local rasters)."""
    geom = ne_gdf.union_all().simplify(tolerance)
    return ee.Geometry(geom.__geo_interface__, None, False)


def run_properties(cfg: Config) -> dict:
    """Identity stamped on every exported asset and checked at download."""
    return {"ch4ne_version": __version__, "ch4ne_start": cfg.start_date, "ch4ne_end": cfg.end_date,
            "ch4ne_label": cfg.run_label, "ch4ne_params": json.dumps(cfg.gee_params(), sort_keys=True)}


def source_collection(cfg: Config, start=None, end=None):
    """L3 orbits in the window that carry every band the pipeline reads."""
    col = (ee.ImageCollection(cfg.collection)
           .filterDate(start or cfg.start_date, end or end_exclusive(cfg))
           .filterBounds(export_region(cfg))
           .filter(ee.Filter.calendarRange(cfg.overpass_start_hour_utc, cfg.overpass_end_hour_utc, "hour")))
    for b in (cfg.band, cfg.uncertainty_band, cfg.aod_band):
        col = col.filter(ee.Filter.listContains("system:band_names", b))
    return col


def collection_versions(cfg: Config):
    """Processor versions in the window and band coverage (the GEE OFFL
    collection mixes processor releases; see docs/AUDITORIA_v3.md, C10)."""
    raw = (ee.ImageCollection(cfg.collection).filterDate(cfg.start_date, end_exclusive(cfg))
           .filterBounds(export_region(cfg)))
    return ee.Dictionary({
        "n_orbits": raw.size(),
        "n_orbits_used": source_collection(cfg).size(),
        "processor_versions": raw.aggregate_histogram("PROCESSOR_VERSION"),
        "first_orbit": ee.Date(raw.aggregate_min("system:time_start")).format("YYYY-MM-dd'T'HH:mm"),
        "last_orbit": ee.Date(raw.aggregate_max("system:time_start")).format("YYYY-MM-dd'T'HH:mm"),
    }).getInfo()


# ------------------------------------------------------------ orbit level
def orbit_collection(cfg: Config):
    proj = grid_projection(cfg)
    t0 = ee.Date(cfg.start_date)

    def prep(img):
        ok = (img.select(cfg.uncertainty_band).lte(cfg.max_uncertainty_ppb)
              .And(img.select(cfg.aod_band).lte(cfg.max_aod_swir)))
        ch4 = img.select(cfg.band).updateMask(ok)
        aod = img.select(cfg.aod_band).updateMask(ch4.mask())
        valid = ch4.mask().gt(0)
        # one reduceResolution, each band weighted by its own mask: mean over
        # the valid 0.01 deg bins of each grid cell (xch4, aod) and share of
        # the cell's bins that hold a valid retrieval (f)
        agg = (ch4.addBands(aod).addBands(valid).rename(["xch4", "aod", "f"])
               .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)
               .reproject(proj))
        agg = agg.select(["xch4", "aod"]).updateMask(agg.select("f").gte(cfg.min_cell_valid_fraction))
        d = img.date()
        return agg.set({
            "month": d.format("YYYY-MM"),
            "month_index": ee.Number(d.get("year")).multiply(12).add(d.get("month")),
            "t_years": d.difference(t0, "year"),
            "janjun": ee.Number(d.get("month")).lte(6),
            "system:time_start": img.get("system:time_start"),
        })

    return source_collection(cfg).map(prep)


# ------------------------------------------------------ monthly background
def _bg_reducer():
    # unweighted: pixel-centre inclusion for BOTH the median and the count;
    # maxRaw above the number of NE cells makes the median exact
    return ee.Reducer.median(maxRaw=60000).combine(ee.Reducer.count(), sharedInputs=True).unweighted()


def month_background_stats(cfg: Config, ne_geom, key, raw=None):
    """NE-land median of the month's mean XCH4 and number of cells.

    Computed from the filtered native orbits evaluated directly on the
    analysis grid: Earth Engine reads the mean-pyramid level closest to
    0.0625 deg, so no per-orbit reduceResolution is needed. B_m is only a
    regional monthly offset (growth, season, processor steps) removed from
    every cell alike; this makes phase A ~10x cheaper than aggregating each
    orbit exactly (measured: 9.7 EECU-hours for the exact version)."""
    raw = raw if raw is not None else filtered_orbits(cfg)
    y, m = map(int, str(key).split("-")) if isinstance(key, str) else (None, None)
    if y is not None:
        start = ee.Date.fromYMD(y, m, 1)
    else:
        start = ee.Date.parse("YYYY-MM", key)
    sub = raw.filterDate(start, start.advance(1, "month"))
    empty = ee.Image.constant(0).rename("xch4").updateMask(0).toFloat()
    mimg = sub.merge(ee.ImageCollection([empty])).mean()   # keeps the band in empty months
    # latitude of the observed cells: a coverage-bias diagnostic (wet months see mostly the south)
    mimg = mimg.addBands(ee.Image.pixelLonLat().select("latitude").updateMask(mimg.mask()))
    st = ee.Dictionary(mimg.reduceRegion(reducer=_bg_reducer(), geometry=ne_geom, crs="EPSG:4326",
                                         crsTransform=crs_transform(cfg), maxPixels=1e9,
                                         tileScale=cfg.tile_scale))
    return st, sub.size()


def filtered_orbits(cfg: Config):
    """Native L3 orbits with the retrieval filters applied, band 'xch4'."""
    def f(img):
        ok = (img.select(cfg.uncertainty_band).lte(cfg.max_uncertainty_ppb)
              .And(img.select(cfg.aod_band).lte(cfg.max_aod_swir)))
        return img.select([cfg.band], ["xch4"]).updateMask(ok)
    return source_collection(cfg).map(f)


def monthly_background(cfg: Config, ne_geom):
    """FeatureCollection(month, bg_ppb, n_cells, n_orbits) + run identity."""
    props = run_properties(cfg)

    raw = filtered_orbits(cfg)

    def one(key):
        key = ee.String(key)
        st, n = month_background_stats(cfg, ne_geom, key, raw)
        # Export.table.toAsset refuses features with null geometry
        return ee.Feature(ee.Geometry.Point([-41.0, -9.0]), {"month": key, "bg_ppb": st.get("xch4_median", NODATA),
                                 "n_cells": st.get("xch4_count", 0), "n_orbits": n,
                                 "lat_median_observed": st.get("latitude_median", NODATA)}).set(props)

    return ee.FeatureCollection(ee.List(month_keys(cfg)).map(one))


def anomaly_collection(cfg: Config, orbits, bg_fc):
    good = bg_fc.filter(ee.Filter.gte("n_cells", cfg.min_bg_cells_per_month))
    months = good.aggregate_array("month")
    bg = ee.Dictionary.fromLists(months, good.aggregate_array("bg_ppb"))

    def add(img):
        b = ee.Number(bg.get(img.get("month")))
        x = img.select("xch4")
        anom = x.subtract(b).rename("anom")
        t = ee.Image.constant(img.get("t_years")).toFloat().updateMask(x.mask()).rename("t")
        jj = ee.Image.constant(img.get("janjun")).toFloat().updateMask(x.mask()).rename("janjun")
        mi = ee.Image.constant(img.get("month_index")).toInt().updateMask(x.mask()).rename("mi")
        return img.addBands([anom, t, jj, mi])

    return orbits.filter(ee.Filter.inList("month", months)).map(add)


# ------------------------------------------------------------ cell level
def covariates(cfg: Config):
    proj = grid_projection(cfg)
    first = ee.Image(ee.ImageCollection("MODIS/061/MCD43A3").first())
    albedo = (ee.ImageCollection("MODIS/061/MCD43A3").filterDate("2022-01-01", "2023-01-01")
              .select("Albedo_BSA_Band7").median().multiply(0.001)
              .setDefaultProjection(first.select("Albedo_BSA_Band7").projection())
              .reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(proj).rename("albedo_b7"))
    dem = ee.Image("CGIAR/SRTM90_V4").select("elevation")
    elev = (dem.reduceResolution(ee.Reducer.mean().combine(ee.Reducer.stdDev(), sharedInputs=True), maxPixels=8192)
            .reproject(proj).rename(["elev_mean", "elev_sd"]))
    occ = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence")
    water = occ.gte(50).unmask(0).rename("w")
    w300 = water.reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(crs="EPSG:4326", scale=300)
    water_frac = w300.reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(proj).rename("water_frac")
    return albedo.addBands(elev).addBands(water_frac).rename(list(COVARIATE_BANDS)).toFloat()


def stats_reducer(cfg: Config):
    """One multi-input reducer, so the prepared orbits are read ONCE.
    Inputs (in order): anom, xch4, t, janjun, aod, mi."""
    anom = (ee.Reducer.mean()
            .combine(ee.Reducer.stdDev(), sharedInputs=True)
            .combine(ee.Reducer.count(), sharedInputs=True)
            .combine(ee.Reducer.percentile([5, 10, 50, 90, 95], maxRaw=cfg.percentile_max_raw),
                     sharedInputs=True))
    return (anom
            .combine(ee.Reducer.mean().forEach(["xch4_mean", "t_mean", "janjun_mean", "aod_mean"]),
                     sharedInputs=False)
            # NonNull: in a multi-input reducer the masked inputs reach countDistinct
            # as null and plain countDistinct counts them (+1 in every cell; checked
            # against the monthly stack on real 2023 data)
            .combine(ee.Reducer.countDistinctNonNull().setOutputs(["n_months"]), sharedInputs=False))


STATS_INPUTS = ["anom", "xch4", "t", "janjun", "aod", "mi"]


def stats_image(cfg: Config, anom):
    """Per-cell statistics of the retrievals (STATS_CORE_BANDS): anomaly mean,
    SD, count and p05/p10/p50/p90/p95 (exact below percentile_max_raw
    retrievals per cell), means of raw XCH4, sampling time, wet-season flag
    and AOD, and the number of distinct months with data. Tail widths come
    from the percentiles. All bands of an orbit share one mask, so the
    multi-input reducer sees exactly the retrievals the anomaly sees."""
    s = anom.select(STATS_INPUTS).reduce(stats_reducer(cfg), cfg.parallel_scale)
    p05, p10, p50, p90, p95 = (s.select(f"p{q}") for q in (5, 10, 50, 90, 95))
    return ee.Image.cat([
        s.select("count"),
        s.select("xch4_mean"),
        s.select("mean"), p50, s.select("stdDev"), p10, p90, p05, p95,
        p90.subtract(p50), p50.subtract(p10), p95.subtract(p50), p50.subtract(p05),
        s.select("t_mean"), s.select("janjun_mean"), s.select("aod_mean"), s.select("n_months"),
    ]).rename(list(STATS_CORE_BANDS)).toFloat()


def monthly_stack(cfg: Config, anom, keys):
    """Monthly mean anomaly (a_YYYY_MM) then monthly count (n_YYYY_MM) for the
    given months, in ONE pass over the orbits: each orbit's anomaly is placed
    in the band of its own month of a fully masked template, so a single
    mean+count reduction fills every month (streaming accumulators, low
    memory). Twelve filtered sub-collection reductions per image were retried
    for lack of memory (measured 2026-09-30)."""
    names = [f"m_{k.replace('-', '_')}" for k in keys]
    template = ee.Image.constant([0] * len(keys)).rename(names).toFloat().updateMask(0)
    band_of = ee.Dictionary(dict(zip(keys, names)))

    def place(img):
        a = img.select(["anom"], ee.List([band_of.get(img.get("month"))]))
        return template.addBands(a, None, True).select(names)

    col = anom.filter(ee.Filter.inList("month", keys)).map(place)
    r = col.reduce(ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True), cfg.parallel_scale)
    means = [r.select(f"{n}_mean").rename(f"a_{n[2:]}") for n in names]
    counts = [r.select(f"{n}_count").rename(f"n_{n[2:]}") for n in names]
    return ee.Image.cat(means + counts).toFloat()


def legacy_image(cfg: Config, ne_geom):
    """Exact v3 composite (for reproducing the published numbers):
    mean of uncertainty-filtered retrievals, reproject(scale=7000), v3 dates."""
    col = (ee.ImageCollection(cfg.collection).filterBounds(ne_geom)
           .filterDate("2020-01-01", "2025-12-31")
           .filter(ee.Filter.calendarRange(cfg.overpass_start_hour_utc, cfg.overpass_end_hour_utc, "hour"))
           .select([cfg.band, cfg.uncertainty_band]))
    return (col.map(lambda i: i.updateMask(i.select(cfg.uncertainty_band).lte(15)))
            .select(cfg.band).mean().reproject(crs="EPSG:4326", scale=cfg.legacy_scale_m)
            .clip(ne_geom).rename("legacy_xch4").toFloat())


# ----------------------------------------------------------------- assets
def asset_root(cfg: Config) -> str:
    suffix = f"_{cfg.run_label}" if cfg.run_label else ""
    return f"projects/{cfg.gee_project}/assets/ch4ne_v4{suffix}"


def asset_id(cfg: Config, name: str) -> str:
    return f"{asset_root(cfg)}/{name}"


def ensure_folder(path: str):
    try:
        ee.data.getAsset(path)
        return
    except ee.EEException:
        pass
    root = path.rsplit("/", 1)[0]
    try:
        ee.data.listAssets({"parent": root})
    except ee.EEException:
        project = root.split("/")[1]
        raise SystemExit(
            f"The Earth Engine asset root '{root}' does not exist yet (usual for a newly registered "
            f"project). Open https://code.earthengine.google.com/?project={project} once, select the "
            "project in the top-right menu and open the 'Assets' tab; this creates the root. Then re-run.")
    ee.data.createAsset({"type": "FOLDER"}, path)


def get_asset(aid):
    try:
        return ee.data.getAsset(aid)
    except ee.EEException:
        return None


def years(cfg: Config):
    return sorted({k[:4] for k in month_keys(cfg)})


def image_jobs(cfg: Config):
    """Image assets of a run. The legacy (v3) composite only makes sense for the main run."""
    names = ["stats"] + [f"monthly_{y}" for y in years(cfg)]
    if cfg.run_label == "" and cfg.is_default_period:
        names.append("legacy")
    return names


def select_jobs(cfg: Config, only=None):
    """Jobs matching `only` (names or prefixes such as 'monthly'); all when None."""
    names = image_jobs(cfg)
    if not only:
        return names
    sel = [n for n in names if n in only or n.split("_")[0] in only]
    unknown = [o for o in only if not any(n == o or n.split("_")[0] == o for n in names)]
    if unknown:
        raise SystemExit(f"Unknown export job(s) {unknown}; available: {names}")
    return sel


def task_description(cfg: Config, name: str) -> str:
    tag = f"_{cfg.run_label}" if cfg.run_label else ""
    return f"ch4ne{tag}_{name}"


def active_tasks(cfg: Config) -> dict:
    """{description: task id} of this project's export tasks that are still queued or running."""
    out = {}
    for op in ee.data.listOperations(f"projects/{cfg.gee_project}"):
        md = op.get("metadata", {})
        if md.get("state") in ACTIVE_STATES:
            out[md.get("description")] = op["name"].rsplit("/", 1)[-1]
    return out


def _reuse_or_start(cfg: Config, name: str, running: dict, overwrite: bool):
    """('running', id) | ('exists', None) | ('start', None) for one asset."""
    desc = task_description(cfg, name)
    if desc in running:
        if overwrite:
            raise SystemExit(f"{desc}: task {running[desc]} is still queued/running. --overwrite would adopt it "
                             "instead of recomputing; cancel it first (https://code.earthengine.google.com/tasks) "
                             "and re-run.")
        print(f"{desc}: a task is already {running[desc]!r} queued/running; waiting for it "
              "(its settings are checked against the configuration at download)")
        return "running", running[desc]
    meta = get_asset(asset_id(cfg, name))
    if meta is not None and not overwrite:
        if name == TABLE_JOB:  # a table carries the identity on its features
            background_table(cfg)
        else:
            _check_identity(cfg, meta, asset_id(cfg, name))
        print(f"{desc}: asset exists with this run's identity; reused")
        return "exists", None
    return "start", None


def start_background_export(cfg: Config, ne_gdf, overwrite=False):
    """Phase A: the monthly NE backgrounds, exported once as a table asset.
    Returns the ids of the tasks to wait for (empty if the asset is reused)."""
    cfg.validate()
    ensure_folder(asset_root(cfg))
    running = active_tasks(cfg)
    if overwrite:
        mine = {task_description(cfg, n) for n in [TABLE_JOB] + image_jobs(cfg)}
        busy = sorted(d for d in running if d in mine)
        if busy:
            raise SystemExit(f"--overwrite: tasks of this run are still queued/running ({busy}); replacing the "
                             "phase-A table now would change the input of running phase-B tasks. Cancel them first "
                             "(https://code.earthengine.google.com/tasks) and re-run.")
    what, tid = _reuse_or_start(cfg, TABLE_JOB, running, overwrite)
    if what == "running":
        return [tid]
    if what == "exists":
        background_table(cfg)  # identity of the table's features
        return []
    bg = monthly_background(cfg, ne_geometry(ne_gdf))
    t = ee.batch.Export.table.toAsset(collection=bg, description=task_description(cfg, TABLE_JOB),
                                      assetId=asset_id(cfg, TABLE_JOB), overwrite=overwrite)
    t.start()
    print("started", t.status()["description"], t.id)
    return [t.id]


def _check_identity(cfg: Config, meta: dict, what: str):
    props = (meta or {}).get("properties", {})
    want = run_properties(cfg)
    got = {k: props.get(k) for k in ("ch4ne_start", "ch4ne_end", "ch4ne_label", "ch4ne_params")}
    exp = {k: want[k] for k in got}
    if got["ch4ne_params"] is not None:
        got["ch4ne_params"] = json.loads(got["ch4ne_params"])
        exp["ch4ne_params"] = json.loads(exp["ch4ne_params"])
    if got != exp:
        raise SystemExit(f"{what} was produced with different settings than the current configuration.\n"
                         f"asset: {got}\nconfig: {exp}\nUse --overwrite to replace it, or another --label.")
    if props.get("ch4ne_version") and props.get("ch4ne_version") != __version__:
        # same settings, other code version: accepted (the version is recorded per asset in
        # PROVENANCE.json), but reported so a change that alters results is not missed
        print(f"NOTE: {what} was produced by ch4ne {props.get('ch4ne_version')}; this code is {__version__}.")


def background_table(cfg: Config):
    aid = asset_id(cfg, TABLE_JOB)
    fc = ee.FeatureCollection(aid)
    first = fc.first().getInfo()
    if first is None:
        raise SystemExit(f"{aid} is empty or missing: run phase A (export) first.")
    _check_identity(cfg, {"properties": first["properties"]}, aid)
    return fc


def _job_image(cfg: Config, name: str, anom, ne_geom):
    """(image, projection kwargs, tile size) of one asset. Tasks run in
    submission order (custom priorities need a commercial project), so the
    order of image_jobs() puts the per-cell statistics first."""
    grid = dict(crs="EPSG:4326", crsTransform=crs_transform(cfg))
    if name == "stats":
        return stats_image(cfg, anom), grid, cfg.stats_tile_px
    if name.startswith("monthly_"):
        keys = [k for k in month_keys(cfg) if k.startswith(name.split("_")[1])]
        return monthly_stack(cfg, anom, keys), grid, cfg.monthly_tile_px
    if name == "legacy":
        return legacy_image(cfg, ne_geom), dict(crs="EPSG:4326", scale=cfg.legacy_scale_m), 256
    raise ValueError(name)


def start_main_exports(cfg: Config, ne_gdf, overwrite=False, only=None):
    """Phase B: per-cell statistics and the yearly monthly stacks,
    reading the phase-A background. Assets that already exist with this run's
    identity, and tasks still running, are reused. Returns task ids to wait for."""
    cfg.validate()
    running = active_tasks(cfg)
    wait_ids, todo = [], []
    for name in select_jobs(cfg, only):
        what, tid = _reuse_or_start(cfg, name, running, overwrite)
        if what == "running":
            wait_ids.append(tid)
        elif what == "start":
            todo.append(name)
    if not todo:
        return wait_ids
    ne_geom = ne_geometry(ne_gdf)
    needs_anom = any(n == "stats" or n.startswith("monthly_") for n in todo)
    anom = anomaly_collection(cfg, orbit_collection(cfg), background_table(cfg)) if needs_anom else None
    props = run_properties(cfg)
    region = export_region(cfg)
    for name in todo:
        img, proj, tile = _job_image(cfg, name, anom, ne_geom)
        # masked pixels become NODATA (pixel downloads would otherwise encode
        # them as 0 ppb); set() comes last so the identity is not dropped
        img = img.unmask(NODATA, False).toFloat().set(props)
        t = ee.batch.Export.image.toAsset(image=img, description=task_description(cfg, name),
                                          assetId=asset_id(cfg, name), region=region, maxPixels=1e10,
                                          overwrite=overwrite, shardSize=tile, **proj)
        t.start()
        print("started", t.status()["description"], t.id, flush=True)
        wait_ids.append(t.id)
    return wait_ids


def task_states(cfg: Config, task_ids):
    """{task_id: (description, state, error)} via the Operations API."""
    out = {}
    for tid in task_ids:
        op = ee.data.getOperation(f"projects/{cfg.gee_project}/operations/{tid}")
        md = op.get("metadata", {})
        out[tid] = (md.get("description", tid), md.get("state", "UNKNOWN"), op.get("error"))
    return out


def optional_descriptions(cfg: Config):
    """Task descriptions whose failure does not stop the run (monthly stacks, v3 composite)."""
    return {task_description(cfg, n) for n in image_jobs(cfg) if n.startswith("monthly_") or n == "legacy"}


def wait(cfg: Config, task_ids, poll_s=60, timeout_h=24, optional=None):
    """Block until every export task finished. Raise if a required task
    failed; failures of optional products (monthly stacks, v3 composite;
    matched by exact task description) are reported and returned, so the core
    analysis can still run. optional=() makes every task required."""
    optional = optional_descriptions(cfg) if optional is None else set(optional)
    if not task_ids:
        return {}
    t_end = time.time() + timeout_h * 3600
    done_states = ("SUCCEEDED", "FAILED", "CANCELLED")
    while True:
        st = task_states(cfg, task_ids)
        print(time.strftime("%H:%M:%S"), {d: s for d, s, _ in st.values()}, flush=True)
        if all(s in done_states for _, s, _ in st.values()):
            break
        if time.time() > t_end:
            raise TimeoutError(st)
        time.sleep(poll_s)
    bad = {d: (s, (e or {}).get("message")) for d, s, e in st.values() if s != "SUCCEEDED"}
    hard = {d: v for d, v in bad.items() if d not in optional}
    if hard:
        raise RuntimeError(f"Export tasks failed: {bad}")
    if bad:
        print("WARNING: optional exports failed:", bad, flush=True)
    return bad


# ----------------------------------------------------------------- checks
def dry_run(cfg: Config, ne_gdf):
    """Quick checks (about a minute): bands, orbit counts, processor versions,
    one prepared orbit and ONE month of background computed directly."""
    cfg.validate()
    first = ee.Image(ee.ImageCollection(cfg.collection).filterDate(cfg.start_date, end_exclusive(cfg))
                     .filterBounds(export_region(cfg)).first())
    print("L3 bands:", first.bandNames().getInfo())
    print("Collection:", json.dumps(collection_versions(cfg), indent=1))
    orbits = orbit_collection(cfg)
    print("Prepared orbit bands:", ee.Image(orbits.first()).bandNames().getInfo())
    print("Orbits with data outside the overpass window (sample days):", overpass_window_leaks(cfg))
    ne_geom = ne_geometry(ne_gdf)
    for key in month_keys(cfg)[:4]:
        st, n = month_background_stats(cfg, ne_geom, key)
        res = ee.Dictionary({"stats": st, "n_orbits": n}).getInfo()
        if res["n_orbits"]:
            print(f"Background {key}:", res)
            break
    else:
        print("No orbits in the first four months of the window.")
    print("Asset folder:", asset_root(cfg))


def overpass_window_leaks(cfg: Config, days=None):
    """Orbits that carry data over the export window but start outside the
    overpass hours (must be 0). Checked on one day per season and year."""
    y0, y1 = int(cfg.start_date[:4]), int(cfg.end_date[:4])
    days = days or [f"{y}-{m:02d}-15" for y in range(y0, y1 + 1) for m in (2, 5, 8, 11)
                    if cfg.start_date <= f"{y}-{m:02d}-15" <= cfg.end_date]
    region = export_region(cfg)
    fc = ee.FeatureCollection([])
    for d in days:
        col = (ee.ImageCollection(cfg.collection).filterDate(d, ee.Date(d).advance(1, "day")).filterBounds(region)
               .filter(ee.Filter.calendarRange(cfg.overpass_start_hour_utc, cfg.overpass_end_hour_utc, "hour").Not()))
        fc = fc.merge(col.map(lambda img: ee.Feature(None, {"n": img.select(cfg.band).reduceRegion(
            ee.Reducer.count(), region, scale=50000, bestEffort=True).get(cfg.band)})))
    return int(fc.filter(ee.Filter.gt("n", 0)).size().getInfo())


# ----------------------------------------------------------------- download
def _fetch_pixels(aid, band_ids, grid):
    data = ee.data.getPixels({"assetId": aid, "fileFormat": "GEO_TIFF", "grid": grid, "bandIds": band_ids})
    if data[:2] == b"PK":  # zipped response
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            data = z.read([n for n in z.namelist() if n.endswith(".tif")][0])
    return data


def _read_tif_bytes(data: bytes):
    import rasterio
    from rasterio.io import MemoryFile
    with MemoryFile(data) as mf, mf.open() as src:
        return src.read(), src.transform, src.width, src.height


def _check_grid(name, transform, width, height, expect):
    ok = (width, height) == (expect["width"], expect["height"]) and         all(abs(a - b) < 1e-9 for a, b in zip(tuple(transform)[:6], expect["transform"]))
    if not ok:
        raise RuntimeError(f"{name}: grid {(width, height, tuple(transform)[:6])} differs from the analysis grid {expect}")


def _write_stack(path: Path, arr, names, transform):
    """Float32 GeoTIFF with NODATA, band descriptions and a band-name sidecar."""
    import rasterio
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[1], width=arr.shape[2], count=arr.shape[0],
                       dtype="float32", crs="EPSG:4326", transform=transform, nodata=NODATA,
                       compress="deflate") as dst:
        dst.write(arr.astype("float32"))
        for i, n in enumerate(names, 1):
            dst.set_band_description(i, n)
    path.with_suffix(".bands.json").write_text(json.dumps(list(names)))
    return path


def _asset_meta(cfg: Config, name: str):
    aid = asset_id(cfg, name)
    meta = get_asset(aid)
    if meta is None:
        return aid, None
    _check_identity(cfg, meta, aid)
    return aid, meta


COVARIATES_FILE = "covariates.tif"


def compute_covariates(cfg: Config, out_dir: Path, block_px=None, threads=4, retries=5):
    """Static covariates on the analysis grid, computed with ee.data.computePixels
    in square blocks of block_px cells (same expression as covariates(), same
    grid, NODATA for masked cells). Written to gee/covariates.tif; reused if it
    already exists for this grid."""
    import numpy as np
    from concurrent.futures import ThreadPoolExecutor

    out = Path(out_dir) / "gee"
    out.mkdir(parents=True, exist_ok=True)
    path = out / COVARIATES_FILE
    w, h = grid_shape(cfg)
    if path.exists():
        arr, tr, ww, hh = _read_tif_bytes(path.read_bytes())
        try:
            _check_grid(COVARIATES_FILE, tr, ww, hh, dict(width=w, height=h, transform=tuple(crs_transform(cfg))))
            if arr.shape[0] == len(COVARIATE_BANDS):
                print(f"{path} exists for this grid; reused")
                return path
        except RuntimeError:
            pass
    b = block_px or cfg.covariates_tile_px
    img = covariates(cfg).unmask(NODATA, False).toFloat()
    res = cfg.grid_res_deg
    blocks = [(r0, c0) for r0 in range(0, h, b) for c0 in range(0, w, b)]
    data = np.full((len(COVARIATE_BANDS), h, w), np.nan, dtype="float32")

    def one(rc):
        r0, c0 = rc
        bh, bw = min(b, h - r0), min(b, w - c0)
        req = {"expression": img, "fileFormat": "NUMPY_NDARRAY",
               "grid": {"dimensions": {"width": bw, "height": bh},
                        "affineTransform": {"scaleX": res, "shearX": 0, "translateX": cfg.grid_west + c0 * res,
                                            "shearY": 0, "scaleY": -res, "translateY": cfg.grid_north - r0 * res},
                        "crsCode": "EPSG:4326"}}
        for k in range(retries):
            try:
                a = ee.data.computePixels(req)
                break
            except ee.EEException as e:
                if k == retries - 1:
                    raise
                print(f"block {rc}: {e}; retrying", flush=True)
                time.sleep(10 * (k + 1))
        for i, n in enumerate(COVARIATE_BANDS):
            data[i, r0:r0 + bh, c0:c0 + bw] = a[n]
        return rc

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=threads) as ex:
        for k, _ in enumerate(ex.map(one, blocks), 1):
            if k % 10 == 0 or k == len(blocks):
                print(f"covariates: {k}/{len(blocks)} blocks ({time.time() - t0:.0f} s)", flush=True)
    if np.isnan(data).any():
        raise RuntimeError("covariates: some blocks were not filled")
    from affine import Affine
    return _write_stack(path, data, COVARIATE_BANDS, Affine(*crs_transform(cfg)))


def write_task_costs(cfg: Config, path: Path):
    """All Earth Engine tasks of this run (by description prefix) with state,
    times and EECU usage, as reported by the Operations API."""
    import pandas as pd
    # every task of this run, including failed or cancelled attempts and earlier job layouts
    # (e.g. the single 'monthly' image, the batch 'covariates' export)
    stems = ("monthly", "stats", "covariates", "legacy", "v4_teste")
    prefix = f"ch4ne_{cfg.run_label}_" if cfg.run_label else "ch4ne_"
    rows = []
    for op in ee.data.listOperations(f"projects/{cfg.gee_project}"):
        md = op.get("metadata", {})
        d = md.get("description", "")
        if d.startswith(prefix) and d[len(prefix):].startswith(stems):
            rows.append(dict(task_id=op["name"].rsplit("/", 1)[-1], description=d, state=md.get("state"),
                             create_time=md.get("createTime"), start_time=md.get("startTime"),
                             end_time=md.get("endTime"), attempt=md.get("attempt"),
                             eecu_hours=(md.get("batchEecuUsageSeconds") or 0) / 3600,
                             error=(op.get("error") or {}).get("message", "")))
    cols = ["task_id", "description", "state", "create_time", "start_time", "end_time", "attempt", "eecu_hours",
            "error"]
    pd.DataFrame(rows, columns=cols).sort_values("create_time").to_csv(path, index=False)
    return path


def download_all(cfg: Config, out_dir: Path):
    """Fetch the finished assets into out_dir/gee and write PROVENANCE.json.

      covariates.tif     = static covariates (computePixels, see compute_covariates)
      stats.tif          = asset 'stats' + covariates.tif (STATS_BANDS order)
      monthly_anom.tif   = a_YYYY_MM bands of every 'monthly_YYYY' asset
      monthly_count.tif  = n_YYYY_MM bands of the same assets
      legacy.tif         = v3 composite (main run only)
      monthly_background.csv

    Stored pixels are read with ee.data.getPixels on the asset's own grid (no
    computation, no time-out; getDownloadURL is avoided because its Python
    client silently drops a camelCase 'crsTransform' key). Values are copied
    unchanged (float32, NODATA for masked cells). The monthly stack is
    optional: if any year is missing it is skipped with a warning (the
    analysis then omits the cluster time series only)."""
    out = Path(out_dir) / "gee"
    out.mkdir(parents=True, exist_ok=True)
    w, h = grid_shape(cfg)
    expect = dict(width=w, height=h, transform=tuple(crs_transform(cfg)))
    files, assets, band_names = {}, {}, {}

    def fetch(name, bands=None):
        aid, meta = _asset_meta(cfg, name)
        if meta is None:
            return None
        ids = [b["id"] for b in meta["bands"]]
        sel = ids if bands is None else [b for b in ids if b in bands or any(b.startswith(x) for x in bands)]
        arr, tr, ww, hh = _read_tif_bytes(_fetch_pixels(aid, sel, meta["bands"][0]["grid"]))
        if arr.shape[0] != len(sel):
            raise RuntimeError(f"{aid}: {arr.shape[0]} bands downloaded, expected {len(sel)}")
        assets[name] = dict(id=meta.get("name", aid), update_time=meta.get("updateTime"),
                            size_bytes=meta.get("sizeBytes"), bands=len(ids),
                            properties={k: v for k, v in meta.get("properties", {}).items()
                                        if k.startswith("ch4ne_")})
        return arr, sel, tr, ww, hh

    # ---- per-cell statistics + covariates -> stats.tif
    import numpy as np
    got = fetch("stats")
    if got is None:
        raise SystemExit(f"{asset_id(cfg, 'stats')} not found: run the export stage first.")
    arr, sel, tr, ww, hh = got
    _check_grid("stats", tr, ww, hh, expect)
    if tuple(sel) != tuple(STATS_CORE_BANDS):
        raise RuntimeError(f"stats: bands {sel} differ from {STATS_CORE_BANDS}")
    files["covariates"] = compute_covariates(cfg, out_dir)
    cov, ctr, cw, ch = _read_tif_bytes(files["covariates"].read_bytes())
    _check_grid(COVARIATES_FILE, ctr, cw, ch, expect)
    band_names["covariates"] = list(COVARIATE_BANDS)
    files["stats"] = _write_stack(out / "stats.tif", np.concatenate([arr, cov]), STATS_BANDS, tr)
    band_names["stats"] = list(STATS_BANDS)

    # ---- monthly stack (one asset per year) -> monthly_anom.tif / monthly_count.tif
    ys = years(cfg)
    have = {y: _asset_meta(cfg, f"monthly_{y}")[1] is not None for y in ys}
    if all(have.values()):
        stacks = {"a_": [], "n_": []}
        names = {"a_": [], "n_": []}
        for y in ys:
            for prefix in ("a_", "n_"):  # two requests per year, far below the 48 MB limit
                arr, sel, tr, ww, hh = fetch(f"monthly_{y}", [prefix])
                _check_grid(f"monthly_{y}", tr, ww, hh, expect)
                stacks[prefix].append(arr)
                names[prefix] += sel
        for part, prefix in (("monthly_anom", "a_"), ("monthly_count", "n_")):
            files[part] = _write_stack(out / f"{part}.tif", np.concatenate(stacks[prefix]), names[prefix], tr)
            band_names[part] = names[prefix]
        if names["a_"] != [f"a_{k.replace('-', '_')}" for k in month_keys(cfg)]:
            raise RuntimeError(f"monthly bands {names['a_']} do not cover the months of the run")
    else:
        print(f"WARNING: monthly assets missing for {[y for y, ok in have.items() if not ok]}; "
              "monthly_anom/monthly_count not written (the cluster time series will be skipped).")
        for part in ("monthly_anom", "monthly_count"):
            for suffix in (".tif", ".bands.json"):
                (out / f"{part}{suffix}").unlink(missing_ok=True)

    # ---- v3 composite (main run only)
    if "legacy" in image_jobs(cfg):
        got = fetch("legacy")
        if got is None:
            print("WARNING: legacy asset missing; the v3 reproduction will be skipped.")
            for suffix in (".tif", ".bands.json"):
                (out / f"legacy{suffix}").unlink(missing_ok=True)
        else:
            arr, sel, tr, ww, hh = got
            files["legacy"] = _write_stack(out / "legacy.tif", arr, sel, tr)
            band_names["legacy"] = sel
    for name in list(files):
        files[f"{name}_bands"] = files[name].with_suffix(".bands.json")

    # ---- monthly background table
    aid = asset_id(cfg, TABLE_JOB)
    fc = background_table(cfg).getInfo()
    import pandas as pd
    bg = pd.DataFrame([f["properties"] for f in fc["features"]])[
        ["month", "bg_ppb", "n_cells", "n_orbits", "lat_median_observed"]]
    bg.sort_values("month").to_csv(out / "monthly_background.csv", index=False)
    files["background"] = out / "monthly_background.csv"
    meta = get_asset(aid)
    assets[TABLE_JOB] = dict(id=meta.get("name", aid), update_time=meta.get("updateTime"),
                             size_bytes=meta.get("sizeBytes"),
                             properties={k: v for k, v in fc["features"][0]["properties"].items()
                                         if k.startswith("ch4ne_")})

    files["task_costs"] = write_task_costs(cfg, out / "gee_tasks.csv")
    write_record(out, files, assets, cfg, extra={
        "collection_versions": collection_versions(cfg),
        "run_properties": run_properties(cfg),   # of this code; each asset's own ones are under assets
        "band_names": band_names,
        "assembly": "stats.tif = asset stats (ee.data.getPixels) + covariates.tif (ee.data.computePixels of "
                    "gee_export.covariates() on the same grid, blocks of covariates_tile_px cells); "
                    "monthly_*.tif = assets monthly_YYYY; values copied unchanged (float32, NODATA=-9999)",
    })
    return files
