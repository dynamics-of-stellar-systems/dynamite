"""chi2 in the ALM loop can come from adelie's residual instead of A @ w.

solve_adelie_alm builds X from A[1:] with unit-L2 column scaling, so
    X[1:] @ beta == A[1:] @ w      (w = beta / col_norm)
and adelie returns state.resid == y - X @ beta. The kinematic part of the
residual is therefore already computed inside the solver, and only row 0 of A
(the total-mass constraint) has to be evaluated separately.

That removes one full pass over A per ALM iteration - the whole matrix, up to
adelie_alm_iters times. Run this file.
"""
import numpy as np
from adelie import solver as _adelie_solver


def demo():
    rng = np.random.default_rng(0)
    n_rows, n_orbs = 400, 25
    A = rng.random((n_rows, n_orbs))
    b = rng.random(n_rows)
    mu = 1e3
    sqrt_mu = np.sqrt(mu)

    A_rest, b_rest = A[1:, :], b[1:]
    X = np.empty((A_rest.shape[0] + 1, n_orbs), order='F')
    X[0, :] = sqrt_mu
    X[1:, :] = A_rest
    col_norm = np.linalg.norm(X, axis=0)
    col_norm[col_norm == 0] = 1.0
    X /= col_norm
    y = np.concatenate([[0.0], b_rest])
    lower, upper = np.zeros(n_orbs), np.full(n_orbs, np.inf)
    weights_arr = np.full(X.shape[0], 1 / X.shape[0])

    lam, state = 0.0, None
    for it in range(6):
        y[0] = sqrt_mu * (1.0 + lam / mu)
        state = _adelie_solver.bvls(X, np.ascontiguousarray(y), lower, upper,
                                    weights=weights_arr, n_threads=1,
                                    tol=1e-10, max_iters=int(2e5),
                                    warm_start=state)
        beta = np.asarray(state.beta).ravel()
        w = (beta / col_norm).astype(np.float64)
        lam -= mu * (w.sum() - 1.0)

        direct = float(np.sum((A @ w - b) ** 2))

        resid = np.asarray(state.resid).ravel()
        row0 = float(A[0] @ w - b[0])
        cheap = row0 * row0 + float(resid[1:] @ resid[1:])

        rel = abs(direct - cheap) / max(abs(direct), 1e-300)
        print(f'  iter {it}: direct={direct:.12e} cheap={cheap:.12e} rel={rel:.2e}')
        assert rel < 1e-9, f'chi2 mismatch at iteration {it}: {rel}'

        # the identity the shortcut rests on
        assert np.allclose(X[1:] @ beta, A_rest @ w)
        assert np.allclose(resid, y - X @ beta)

    print("ALM chi2 from adelie residual matches A @ w, OK")


if __name__ == "__main__":
    demo()
