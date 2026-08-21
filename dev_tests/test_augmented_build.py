"""Bit-identity of the extracted augmented-matrix builder vs the inline code
it replaces in solve_adelie_alm (weight_solvers.py).

Run from the repo root:
    PYTHONPATH=. python dev_tests/test_augmented_build.py
"""

import numpy as np

from dynamite.weight_solvers import NNLS


def _reference(A_rest, b_rest, sqrt_mu, dtype):
    """Transcription of the pre-extraction inline construction."""
    X = np.empty((A_rest.shape[0] + 1, A_rest.shape[1]), dtype=dtype, order="F")
    X[0, :] = sqrt_mu
    X[1:, :] = A_rest
    col_norm = np.linalg.norm(X, axis=0)
    col_norm[col_norm == 0] = 1.0
    X /= col_norm
    y = np.concatenate([[0.0], b_rest]).astype(dtype)
    return X, col_norm, y


def test_bitwise():
    rng = np.random.default_rng(42)
    for dtype in (np.float64, np.float32):
        A_rest = rng.standard_normal((37, 11)).astype(dtype)
        A_rest[:, 3] = 0.0  # exercise zero-column guard
        b_rest = rng.standard_normal(37).astype(dtype)
        sqrt_mu = np.sqrt(1e7)
        got = NNLS._build_augmented_X(A_rest, b_rest, sqrt_mu, dtype)
        ref = _reference(A_rest, b_rest, sqrt_mu, dtype)
        for g, r in zip(got, ref):
            assert g.dtype == r.dtype and g.shape == r.shape
            assert np.array_equal(g, r), dtype


if __name__ == "__main__":
    test_bitwise()
    print("test_augmented_build OK")
