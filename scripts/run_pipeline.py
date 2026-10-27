"""
Command-line driver for the v4 pipeline (run from the repository root).

    python scripts/run_pipeline.py check      # Earth Engine login + quick checks (about a minute)
    python scripts/run_pipeline.py export     # phase A (monthly background) -> wait -> phase B
                                              # (stats, monthly_YYYY, legacy); covariates are
                                              # computed at download (computePixels)
    python scripts/run_pipeline.py wait       # poll the phase-B tasks until they finish
    python scripts/run_pipeline.py download   # fetch GeoTIFF/CSV + write PROVENANCE.json
    python scripts/run_pipeline.py analyze    # tables, figures, numbers, checks
    python scripts/run_pipeline.py all        # export -> wait -> download -> analyze

    # robustness runs on the two processor periods (labels are mandatory for non-default dates);
    # --min-obs keeps the main rule's density (50 retrievals in 72 months) for the shorter periods
    python scripts/run_pipeline.py all --start 2020-01-01 --end 2022-07-16 --label P1 --min-obs 22
    python scripts/run_pipeline.py all --start 2022-07-17 --end 2025-12-31 --label P2 --min-obs 29
    python scripts/run_pipeline.py compare --compare P1 P2

    # another Google Cloud project with Earth Engine enabled
    python scripts/run_pipeline.py check --project my-project-id

    # export only some assets (names or prefixes); finished assets are always reused
    python scripts/run_pipeline.py export --jobs stats monthly

Exports are resumable: an asset that already exists with this run's identity
is reused and a task still running is waited for, so re-running a stage never
repeats finished work (use --overwrite to recompute).

Outputs go to outputs_v4[_LABEL]/ (rasters in .../gee/). The export stage
writes run_config.json; later stages reuse those settings and stop if the
command line asks for something different.
"""
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ch4ne.config import Config  # noqa: E402


def build_config(args) -> Config:
    cfg = Config(base_dir=ROOT, run_label=args.label)
    if args.project:
        cfg.gee_project = args.project
    if args.start:
        cfg.start_date = args.start
    if args.end:
        cfg.end_date = args.end
    if args.n_boot:
        cfg.n_boot = args.n_boot
    if args.min_obs:
        cfg.min_obs = args.min_obs
    return cfg.validate()


def manifest_path(cfg):
    return cfg.out_dir / "run_config.json"


def write_manifest(cfg, phase_a, phase_b):
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    p = manifest_path(cfg)
    old = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    history = old.get("task_history", [])
    for tid in list(old.get("tasks_phase_a", [])) + list(old.get("tasks_phase_b", [])) + phase_a + phase_b:
        if tid not in history:
            history.append(tid)
    p.write_text(json.dumps(dict(
        created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        gee_params=cfg.gee_params(), gee_project=cfg.gee_project, config=cfg.as_dict(),
        tasks_phase_a=phase_a, tasks_phase_b=phase_b, task_history=history), indent=2), encoding="utf-8")


def check_manifest(cfg):
    """Later stages must use the settings of the export that produced the assets."""
    p = manifest_path(cfg)
    if not p.exists():
        return None
    m = json.loads(p.read_text(encoding="utf-8"))
    if m["gee_params"] != json.loads(json.dumps(cfg.gee_params())):
        raise SystemExit(f"{p} was written by an export with different settings:\n{m['gee_params']}\n"
                         f"current:\n{cfg.gee_params()}\nUse the same --start/--end/--label, or re-export.")
    if m.get("gee_project") and m["gee_project"] != cfg.gee_project:
        raise SystemExit(f"{p} was exported with project {m['gee_project']}; pass --project {m['gee_project']}.")
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["check", "export", "wait", "download", "analyze", "all", "compare"])
    ap.add_argument("--label", default="", help="run label; separates outputs and GEE assets (e.g. P1, P2)")
    ap.add_argument("--start", default=None, help="first date YYYY-MM-DD (default 2020-01-01)")
    ap.add_argument("--end", default=None, help="last date YYYY-MM-DD, inclusive (default 2025-12-31)")
    ap.add_argument("--project", default=None, help="Google Cloud project with Earth Engine "
                    "(default: CH4NE_GEE_PROJECT or ch4ne/config.py)")
    ap.add_argument("--overwrite", action="store_true", help="recompute and replace existing assets of this label")
    ap.add_argument("--jobs", nargs="+", default=None, metavar="JOB",
                    help="export only these assets: stats, monthly (all years) or monthly_YYYY, legacy")
    ap.add_argument("--no-sensitivity", action="store_true", help="skip the sensitivity grid (faster)")
    ap.add_argument("--n-boot", type=int, default=None, help="bootstrap replicates (default from config)")
    ap.add_argument("--min-obs", type=int, default=None,
                    help="minimum retrievals per cell (analysis only; default from config). For a shorter period "
                         "keep the same density per month, e.g. 50 * n_months / 72")
    ap.add_argument("--compare", nargs=2, metavar=("LABEL_A", "LABEL_B"), help="labels for the compare stage")
    args = ap.parse_args()

    if args.stage == "compare":
        if not args.compare:
            raise SystemExit("compare needs --compare LABEL_A LABEL_B")
        from ch4ne.analysis import compare_runs
        a, b = args.compare
        oa, ob = Config(base_dir=ROOT, run_label=a).out_dir, Config(base_dir=ROOT, run_label=b).out_dir
        dest = Config(base_dir=ROOT).out_dir / "tables" / f"TableS12_clusters_{a}_vs_{b}.csv"
        df = compare_runs(oa, ob, dest)
        print(df.to_string(index=False))
        print("Written", dest)
        return

    cfg = build_config(args)
    gee_dir = cfg.out_dir / "gee"

    if args.stage in ("check", "export", "wait", "download", "all"):
        from ch4ne import gee_export as G
        from ch4ne.boundaries import load_states, ne_states
        G.initialize(cfg)
        ne = ne_states(load_states(cfg), cfg)

    if args.stage == "check":
        G.dry_run(cfg, ne)
        return

    if args.stage in ("export", "all"):
        G.select_jobs(cfg, [j for j in (args.jobs or []) if j != "monthly_background"] or None)   # reject unknown --jobs early
        # phase A is replaced only by a full --overwrite (no --jobs) or when explicitly selected;
        # otherwise selected phase-B assets would be recomputed against a table the others did not use
        redo_a = args.overwrite and (not args.jobs or "monthly_background" in args.jobs)
        ta = G.start_background_export(cfg, ne, overwrite=redo_a)
        write_manifest(cfg, ta, [])
        G.wait(cfg, ta, optional=())              # phase B reads the phase-A asset (never optional)
        b_jobs = [j for j in (args.jobs or []) if j != "monthly_background"]
        tb = ([] if (args.jobs and not b_jobs)
              else G.start_main_exports(cfg, ne, overwrite=args.overwrite, only=b_jobs or None))
        write_manifest(cfg, ta, tb)
    else:
        check_manifest(cfg)

    if args.stage in ("wait", "all"):
        m = check_manifest(cfg)
        if m is None:
            raise SystemExit(f"{manifest_path(cfg)} not found: run the 'export' stage first.")
        ids = list(m.get("tasks_phase_b", []))
        # tasks of earlier exports of this run (e.g. other --jobs) that are still active
        active = set(G.active_tasks(cfg).values())
        ids += [t for t in m.get("task_history", []) if t in active and t not in ids]
        G.wait(cfg, ids)

    if args.stage in ("download", "all"):
        files = G.download_all(cfg, cfg.out_dir)
        print("Downloaded:", {k: str(v) for k, v in files.items()})

    if args.stage in ("analyze", "all"):
        from ch4ne.analysis import run
        res = run(cfg, gee_dir, cfg.out_dir, sensitivity=not args.no_sensitivity)
        print(json.dumps(res["checks"], indent=2, default=str))
        print("Outputs in", res["out"])


if __name__ == "__main__":
    main()
