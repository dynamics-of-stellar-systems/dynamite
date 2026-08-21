"""Real-library proof that the fused adelie constructor is bitwise identical
to the two-step construction it replaces.

Two invocations, each peaking at ~1x matrix + orblib histograms (never 2x),
so they are safe to run one after the other even at production scale:

    ... _real_fused_check.py CONFIG --mode reference --out ref.json
    ... _real_fused_check.py CONFIG --mode fused     --out fused.json

Any diff of the two JSON files (excluding 'timing_s') must be empty: X is
digested as sha256 of contiguous row slabs, so equality proves bit-identity;
col_norm/y/row0_vec are compared through the hex of their exact bytes.

Run from PM_grid (config paths resolve against cwd):
    OPENBLAS_CORETYPE=Haswell PYTHONPATH=<repo> \
    ENV/bin/python <repo>/dev_tests/_real_fused_check.py \
      NGC5139_config_adelie_xeast_profile.yaml --mode reference --out fused_ref.json
"""

import argparse
import hashlib
import json
import time

import numpy as np

import dynamite as dyn


def sha(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def digest_X(X, slab):
    return [sha(X[i : i + slab]) for i in range(0, X.shape[0], slab)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--mode", required=True, choices=["reference", "fused"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--slab", type=int, default=2048)
    args = ap.parse_args()

    c = dyn.config_reader.Configuration(args.config, reset_logging=True)
    model = c.all_models.get_model_from_row(0)
    solver = dyn.weight_solvers.NNLS(config=c, model=model)
    orblib = model.get_orblib()
    t0 = time.perf_counter()
    orblib.read_vel_histograms()
    print(f"orblib read {time.perf_counter() - t0:.1f}s", flush=True)

    NNLS = dyn.weight_solvers.NNLS
    if args.mode == "reference":
        t0 = time.perf_counter()
        A, b = solver.construct_nnls_matrix_and_rhs(orblib)
        X, col_norm, y = NNLS._build_augmented_X(
            A[1:], b[1:], np.sqrt(solver.adelie_mu), A.dtype
        )
        row0_vec = A[0]
        b0 = float(b[0])
        del A, b
        timing = time.perf_counter() - t0
    else:
        t0 = time.perf_counter()
        problem = solver.construct_adelie_matrix_and_rhs(orblib)
        X = problem.X
        col_norm = problem.col_norm
        y = problem.y
        row0_vec = problem.row0_vec
        b0 = problem.b0
        timing = time.perf_counter() - t0

    out = {
        "mode": args.mode,
        "timing_s": timing,
        "shape": list(X.shape),
        "dtype": str(X.dtype),
        "b0": b0,
        "col_norm_sha": sha(col_norm),
        "y_sha": sha(y),
        "row0_sha": sha(row0_vec),
        "x_slab_rows": args.slab,
        "x_slabs": digest_X(X, args.slab),
    }
    del X
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(
        f"wrote {args.out}: {len(out['x_slabs'])} slabs, "
        f"X {out['shape']} {out['dtype']}, build {timing:.1f}s"
    )


if __name__ == "__main__":
    main()
