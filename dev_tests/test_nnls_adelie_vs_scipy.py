#!/usr/bin/env python3
"""Run the NGC6278 weight-solving test with BOTH nnls_solvers and compare.

This is the weight-solver half of test_nnls.py. It reuses the orbit libraries
already in NGC6278_output rather than re-integrating (orblib_f_new_mirror is
not built on this machine), and checks each solver against the same reference
chi2 values test_nnls.py uses: data/chi2_compare_ml_1086.dat.

For each model it reports chi2, kinchi2, sum(w) and the KKT violation, so a
solver that silently returns a non-optimal point is visible.

The models' orbit_weights.ecsv files are backed up and restored, so the tree is
left exactly as found.

Usage:
    python test_nnls_adelie_vs_scipy.py [--config user_test_config_ml.yaml]
"""

import argparse
import os
import shutil
import sys
import time

import numpy as np
from astropy import table

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import dynamite as dyn
from dynamite.weight_solvers import NNLS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="user_test_config_ml.yaml")
    ap.add_argument("--solvers", default="scipy,adelie")
    ap.add_argument("--compare",
                    default="data/chi2_compare_ml_1086.dat")
    args = ap.parse_args()

    print(f"dynamite {dyn.__version__} from {dyn.__path__[0]}", flush=True)

    ref = None
    if os.path.isfile(args.compare):
        ref = table.Table.read(args.compare, format="ascii")
        print(f"reference: {args.compare}", flush=True)

    results = {}
    for solver in [s.strip() for s in args.solvers.split(",")]:
        c = dyn.config_reader.Configuration(
            args.config, reset_logging=True,
            user_logfile=f"test_nnls_{solver}", reset_existing_output=False)
        c.settings.weight_solver_settings["nnls_solver"] = solver
        rows = []
        for i in range(len(c.all_models.table)):
            mod = c.all_models.get_model_from_row(i)
            wfile = os.path.join(mod.directory, dyn.constants.weight_file)
            if not os.path.isdir(mod.directory):
                continue
            backup = wfile + ".bak_solvertest"
            had = os.path.isfile(wfile)
            if had:
                shutil.copy2(wfile, backup)
            try:
                orblib = mod.get_orblib()
                ws = NNLS(config=c, model=mod)
                t0 = time.perf_counter()
                # ignore_existing_weights: force a real solve, not a file read
                w, chi2, kinchi2, kinmap = ws.solve(
                    orblib, ignore_existing_weights=True)
                dt = time.perf_counter() - t0
                orblib.read_vel_histograms()
                A, b = ws.construct_nnls_matrix_and_rhs(orblib)
                kkt, kkt_raw = NNLS.kkt_violation(A, b, w)
                rows.append(dict(i=i, ml=float(c.all_models.table["ml"][i]),
                                 chi2=chi2, kinchi2=kinchi2, t=dt,
                                 sw=float(np.sum(w)),
                                 act=int(np.sum(w > 0)), kkt=kkt,
                                 kkt_raw=kkt_raw))
                print(f"  [{solver}] model {i} ml={rows[-1]['ml']:.2f}: "
                      f"chi2={chi2:.2f} kinchi2={kinchi2:.2f} "
                      f"sum_w={rows[-1]['sw']:.8f} act={rows[-1]['act']} "
                      f"KKT={kkt:.2e} raw={kkt_raw:.2e} ({dt:.1f}s)", flush=True)
            except Exception as e:
                print(f"  [{solver}] model {i}: FAILED "
                      f"{type(e).__name__}: {e}", flush=True)
            finally:
                if had:
                    shutil.move(backup, wfile)
                elif os.path.isfile(wfile):
                    os.remove(wfile)
        results[solver] = rows

    print("\n" + "=" * 104)
    print(f"{'model':>5s} {'ml':>6s} " +
          "".join(f"{s + ' chi2':>15s}{s + ' kinchi2':>16s}"
                  for s in results) + f"{'ref chi2':>12s}")
    print("-" * 104)
    solvers = list(results)
    n = max(len(v) for v in results.values()) if results else 0
    for k in range(n):
        line = ""
        ml = None
        for s in solvers:
            r = results[s][k] if k < len(results[s]) else None
            if r:
                ml = r["ml"]
                line += f"{r['chi2']:15.2f}{r['kinchi2']:16.2f}"
            else:
                line += f"{'--':>15s}{'--':>16s}"
        rv = ""
        if ref is not None and k < len(ref):
            rv = f"{float(ref['chi2'][k]):12.2f}"
        print(f"{k:>5d} {ml if ml else 0:6.2f} {line}{rv}")
    print("=" * 104)

    if len(solvers) == 2 and all(results[s] for s in solvers):
        a, bb = results[solvers[0]], results[solvers[1]]
        print(f"\n{solvers[1]} / {solvers[0]} ratios:")
        for x, y in zip(a, bb):
            print(f"  ml={x['ml']:.2f}  chi2 {y['chi2']/x['chi2']:.6f}   "
                  f"kinchi2 {y['kinchi2']/x['kinchi2']:.6f}   "
                  f"KKT {x['kkt']:.1e} -> {y['kkt']:.1e}   "
                  f"speed {x['t']/max(y['t'],1e-9):.2f}x")


if __name__ == "__main__":
    main()
