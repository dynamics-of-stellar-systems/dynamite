"""End-to-end check of the ALM chi2 shortcut against a real NNLS matrix.

solve_adelie_alm returns the best iterate and logs its chi2. That chi2 is now
derived from adelie's residual instead of A @ w, so recomputing it directly
from the returned weights is an independent check of the identity.

Uses a row-subsampled slice of a real A (row 0, the total-mass constraint, is
always kept) so it finishes on a laptop.

Run from ~/research/omegacen/dynamite_dataprep with the `main` conda python:
    python .../_real_alm_chi2_check.py /tmp/old.npz [n_rows]
"""

import logging
import sys

import numpy as np

import dynamite as dyn

CONFIG = "NGC5139_config_veldist_combined_bigger.yaml"


def main(npz, n_rows=40000):
    d = np.load(npz)
    A_full, b_full = d["A"], d["b"]
    rng = np.random.default_rng(0)
    keep = np.concatenate(
        [[0], 1 + rng.choice(A_full.shape[0] - 1, n_rows - 1, replace=False)]
    )
    keep.sort()
    A = np.asfortranarray(A_full[keep])
    b = np.ascontiguousarray(b_full[keep])
    del A_full, b_full
    print(f"A {A.shape} {A.dtype}, F-contiguous={A.flags.f_contiguous}", flush=True)

    c = dyn.config_reader.Configuration(CONFIG, reset_logging=True)
    model = c.all_models.get_model_from_row(0)
    solver = dyn.weight_solvers.NNLS(config=c, model=model)
    solver.nnls_solver = "adelie"

    # capture the chi2 the solver reports, which is the residual-derived one
    captured = {}

    class _Grab(logging.Handler):
        def emit(self, record):
            msg = record.getMessage()
            if "chi2=" in msg:
                captured["chi2"] = float(msg.split("chi2=")[1].split(",")[0])

    log = logging.getLogger("dynamite.weight_solvers.NNLS")
    log.setLevel(logging.INFO)
    log.addHandler(_Grab())

    problem = _problem_from_A_b(solver, A, b)
    w, resid_full = solver.solve_adelie_alm(problem)

    # recompute chi2 for the returned iterate the expensive, direct way
    direct = float(np.sum((A.astype(np.float64) @ w - b.astype(np.float64)) ** 2))
    reported = captured.get("chi2")
    # and the fused chi2_vector path solve() now uses
    row0_resid = float(problem.row0_vec @ w) - problem.b0
    print(f"\nchi2 reported by solver (from resid) : {reported!r}")
    print(f"chi2 recomputed directly as A@w - b  : {direct:.10f}")
    if reported is not None:
        rel = abs(direct - reported) / abs(direct)
        print(
            f"relative difference                  : {rel:.2e}"
            f"   (solver logs 4 dp, so ~1e-9 is exact agreement)"
        )
    print(f"sum(w) = {w.sum():.10f}")
    return direct, reported


def _problem_from_A_b(solver, A, b):
    """Build an AdelieProblem from a classic (A, b) pair via the shared
    finishing steps — the same thing construct_adelie_matrix_and_rhs does
    without ever materializing A."""
    sqrt_mu = np.sqrt(solver.adelie_mu)
    dtype = A.dtype
    X, col_norm, y = dyn.weight_solvers.NNLS._build_augmented_X(
        A[1:], b[1:], sqrt_mu, dtype
    )
    return dyn.weight_solvers.AdelieProblem(
        X=X, col_norm=col_norm, y=y, row0_vec=A[0], b0=float(b[0])
    )


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 40000)
