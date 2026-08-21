"""Profile the RSS timeline of a production-shape NNLS.solve(), phase by phase.

Drives the real pipeline exactly like solve() does (orblib read -> matrix ->
ALM -> final chi2), wrapping each stage so a background RssSampler tags the
timeline. Also counts ALM iterations via the bvls entry point.

The --config yaml must point at a SCRATCH copy of the baseline output tree so
solve()'s weight-file write cannot touch recorded evidence. --alm-iters caps
the multiplier loop for fast memory observation only (per-iteration memory is
flat); never use capped runs for numerical comparisons.

Run (from the dynamite repo root):
    OPENBLAS_CORETYPE=Haswell OMP_NUM_THREADS=24 \
    PYTHONPATH=$PWD \
    ENV/bin/python dev_tests/_real_solve_rss_profile.py \
      /nexus/.../PM_grid/NGC5139_config_adelie_xeast_profile.yaml \
      --tag base_f64_cap30 --alm-iters 30
"""

import argparse
import json
import resource
import sys
import time

import numpy as np

import dynamite as dyn
import dynamite.weight_solvers as ws
from _rss_probe import RssSampler

PM_GRID = "/nexus/posix0/MIA-astro-env/nneum/pesmith/PM_grid"


def rss_gib():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == "darwin" else r / 2**20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--alm-iters", type=int, default=None)
    ap.add_argument("--dtype", default=None, choices=["float32", "float64"])
    args = ap.parse_args()

    csv_path = f"{PM_GRID}/rss_{args.tag}.csv"

    c = dyn.config_reader.Configuration(args.config, reset_logging=True)
    model = c.all_models.get_model_from_row(0)
    solver = dyn.weight_solvers.NNLS(config=c, model=model)
    solver.nnls_solver = "adelie"
    if args.dtype == "float32":
        solver.nnls_dtype = np.float32
    if args.alm_iters is not None:
        solver.adelie_alm_iters = args.alm_iters

    orblib = model.get_orblib()
    phases = {}
    cur = {"name": "startup"}
    samp = RssSampler(csv_path)

    def mark(name):
        phases[name] = {"t": time.perf_counter(), "peak_rss_gib": rss_gib()}
        print(
            f"[phase] {name}: t={phases[name]['t']:.1f}s "
            f"peak_rss={phases[name]['peak_rss_gib']:.1f}GiB",
            flush=True,
        )
        cur["name"] = name
        samp.phase = name

    # wrap instance methods so the real call graph is what gets attributed
    def wrap(obj, name, label):
        orig = getattr(obj, name)

        def wrapper(*a, **k):
            mark(f"{label}:enter")
            try:
                return orig(*a, **k)
            finally:
                mark(f"{label}:exit")

        setattr(obj, name, wrapper)

    wrap(orblib, "read_vel_histograms", "orblib_read")
    ctor = (
        "construct_adelie_matrix_and_rhs"
        if hasattr(type(solver), "construct_adelie_matrix_and_rhs")
        else "construct_nnls_matrix_and_rhs"
    )
    wrap(solver, ctor, "matrix_build")
    wrap(solver, "solve_adelie_alm", "alm_solve")
    kkt_name = (
        "kkt_violation_augmented"
        if hasattr(type(solver), "kkt_violation_augmented")
        else "kkt_violation"
    )
    wrap(solver, kkt_name, "kkt")

    # count ALM iterations without touching library code
    alm_calls = {"n": 0}
    orig_bvls = ws._adelie_solver.bvls

    def counting_bvls(*a, **k):
        alm_calls["n"] += 1
        if alm_calls["n"] % 10 == 1:
            mark(f"alm_iter_{alm_calls['n']}")
        return orig_bvls(*a, **k)

    ws._adelie_solver.bvls = counting_bvls

    t0 = time.perf_counter()
    with samp:
        mark("startup")
        weights, chi2_tot, chi2_kin, chi2_kinmap = solver.solve(
            orblib, ignore_existing_weights=True
        )
        samp.phase = "done"
    wall = time.perf_counter() - t0
    ws._adelie_solver.bvls = orig_bvls

    summary = {
        "tag": args.tag,
        "wall_s": wall,
        "peak_rss_gib": rss_gib(),
        "alm_iters_run": alm_calls["n"],
        "chi2_tot": None if np.isnan(chi2_tot) else float(chi2_tot),
        "chi2_kin": None if np.isnan(chi2_kin) else float(chi2_kin),
        "phases": phases,
    }
    with open(csv_path.replace(".csv", "_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    np.save(f"{PM_GRID}/weights_{args.tag}.npy", weights)


if __name__ == "__main__":
    main()
