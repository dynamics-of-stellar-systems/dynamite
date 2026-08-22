"""Surrogate chi2-vector and KKT identities vs direct computation.

The surrogates are algebraically identical to the A-based forms but re-order
the floating-point operations, so they are tested to rtol 1e-10, NOT bitwise.
The scaled KKT value must stay inside [0, 1].

Run from the repo root:
    PYTHONPATH=. python dev_tests/test_surrogate_chi2_kkt.py
"""

import numpy as np
from scipy.optimize import nnls

from dynamite.weight_solvers import NNLS, chi2_vector_from_residuals


def _problem(rng, n_rows, n_orbs, dtype=np.float64):
    A = rng.standard_normal((n_rows, n_orbs)).astype(dtype)
    A[0, :] = 1e8  # mimic the total-mass row scale
    w_true = rng.random(n_orbs) ** 3  # mostly-near-zero weights
    w_true[:3] = 0.0  # force some exact zeros
    b = A @ w_true.astype(dtype)
    return A, b


def test_kkt_matches_stock():
    rng = np.random.default_rng(7)
    A, b = _problem(rng, 300, 40)
    mu = 1e7
    X, col_norm, _y = NNLS._build_augmented_X(A[1:], b[1:], np.sqrt(mu), np.float64)
    row0_vec = A[0]  # the true total-mass row
    w, _ = nnls(A, b)
    resid_full = A @ w - b  # plain residual aligned to A rows
    got = NNLS.kkt_violation_augmented(row0_vec, b[0], X, col_norm, resid_full, w, mu)
    ref = NNLS.kkt_violation(A, b, w)
    assert 0.0 <= got[0] <= 1.0, got
    assert np.isclose(got[0], ref[0], rtol=1e-10), (got, ref)
    assert np.isclose(got[1], ref[1], rtol=1e-8), (got, ref)


def test_exact_fit_returns_zero_scaled():
    rng = np.random.default_rng(11)
    n_orbs, n_mass = 12, 30
    A = rng.standard_normal((n_mass, n_orbs))
    w = np.zeros(n_orbs)
    w[:5] = rng.random(5) + 0.1
    b = A @ w  # exactly attainable
    mu = 1e7
    X, col_norm, _y = NNLS._build_augmented_X(A[1:], b[1:], np.sqrt(mu), np.float64)
    resid_full = A @ w - b  # ~0 everywhere -> ||r|| ~ 0 guard
    scaled, raw = NNLS.kkt_violation_augmented(
        A[0], b[0], X, col_norm, resid_full, w, mu
    )
    assert raw <= 1e-8 and scaled <= 1e-8, (scaled, raw)


def test_chi2_vector_identity():
    rng = np.random.default_rng(13)
    A, b = _problem(rng, 120, 25)
    w, _ = nnls(A, b)
    resid_full = A @ w - b
    got = chi2_vector_from_residuals(resid_full, float(resid_full[0]) ** 2)
    ref = (A @ w - b) ** 2
    assert np.array_equal(got, ref)  # pure reshape/square: bitwise
    assert got.shape == ref.shape


if __name__ == "__main__":
    test_kkt_matches_stock()
    test_exact_fit_returns_zero_scaled()
    test_chi2_vector_identity()
    print("test_surrogate_chi2_kkt OK")
