"""The fused constructor's zero-error guard: row-restricted == full mask.

construct_adelie_matrix_and_rhs used to find nonzero matrix coefficients on
zero-error rows with a full (n_rows-1, n_orbs) boolean mask, which is ~16 GiB
at production scale (and ~31 GiB with the `&` temporary). It now restricts to
the zero-error rows first. This checks the two agree - which rows/columns are
reported, and which entries get zeroed - including the edge cases that decide
whether the guard fires at all.

Run: python test_zero_error_guard.py
"""

import numpy as np


def old_guard(X, econ_body):
    """The full-mask version, transcribed from before the change."""
    X = X.copy()
    nz = (X[1:, :] != 0) & (econ_body == 0)[:, None]
    fired = bool(np.any(nz))
    reported = None
    if fired:
        rr, oo = np.nonzero(nz)
        reported = (rr + 1, oo, X[1:, :][nz].copy())
        X[1:, :][nz] = 0
    return X, fired, reported


def new_guard(X, econ_body):
    """The row-restricted version now in weight_solvers.py."""
    X = X.copy()
    fired = False
    reported = None
    bad = np.nonzero(econ_body == 0)[0]
    if bad.size:
        blk = X[1:, :][bad]
        nz_rows, nz_cols = np.nonzero(blk)
        if nz_rows.size:
            fired = True
            rr = bad[nz_rows]
            reported = (rr + 1, nz_cols, blk[nz_rows, nz_cols].copy())
            X[1 + rr, nz_cols] = 0
    return X, fired, reported


def check(name, X, econ_body):
    Xo, fired_o, rep_o = old_guard(X, econ_body)
    Xn, fired_n, rep_n = new_guard(X, econ_body)
    assert fired_o == fired_n, f"{name}: guard fired {fired_o} vs {fired_n}"
    assert np.array_equal(Xo, Xn), f"{name}: zeroed matrices differ"
    if rep_o is None:
        assert rep_n is None, f"{name}: only one version reported"
    else:
        # both must name the same [constraint no, orbit no] pairs and values.
        # np.nonzero on the full mask is row-major over all rows; the
        # restricted one is row-major over the bad rows only - same order,
        # since bad rows are ascending. Sort anyway so the check is about the
        # content, not the traversal.
        for i, what in enumerate(("rows", "cols", "values")):
            a = np.sort(np.asarray(rep_o[i]))
            b = np.sort(np.asarray(rep_n[i]))
            assert np.array_equal(a, b), f"{name}: reported {what} differ"
    print(f"  {name}: fired={fired_o}, OK")


def main():
    rng = np.random.default_rng(20260823)
    n_rows, n_orbs = 40, 7

    # 1. the normal case: no zero-error rows at all, guard must not fire
    X = rng.normal(size=(n_rows, n_orbs))
    econ = np.abs(rng.normal(size=n_rows - 1)) + 0.5
    check("no zero-error rows", X, econ)

    # 2. zero-error rows, but the matrix is zero there - guard must not fire.
    #    This is the case the guard exists to let through.
    X = rng.normal(size=(n_rows, n_orbs))
    econ = np.abs(rng.normal(size=n_rows - 1)) + 0.5
    econ[[3, 11, 29]] = 0.0
    X[1 + np.array([3, 11, 29]), :] = 0.0
    check("zero-error rows, zero coeffs", X, econ)

    # 3. zero-error rows with nonzero coefficients - guard must fire
    X = rng.normal(size=(n_rows, n_orbs))
    econ = np.abs(rng.normal(size=n_rows - 1)) + 0.5
    econ[[3, 11, 29]] = 0.0
    check("zero-error rows, nonzero coeffs", X, econ)

    # 4. partially zero row: only some columns offend
    X = rng.normal(size=(n_rows, n_orbs))
    econ = np.abs(rng.normal(size=n_rows - 1)) + 0.5
    econ[7] = 0.0
    X[8, :] = 0.0
    X[8, [0, 4]] = [1.5, -2.5]
    check("partially zero offending row", X, econ)

    # 5. every row has zero error
    X = rng.normal(size=(n_rows, n_orbs))
    check("all rows zero-error", X, np.zeros(n_rows - 1))

    # 6. first and last body rows offend (boundary of the 1+ offset)
    X = rng.normal(size=(n_rows, n_orbs))
    econ = np.abs(rng.normal(size=n_rows - 1)) + 0.5
    econ[[0, n_rows - 2]] = 0.0
    check("first and last body rows", X, econ)

    # 7. row 0 of X is the ALM penalty row and has no econ semantics: a
    #    nonzero there must never be touched by either version
    X = rng.normal(size=(n_rows, n_orbs))
    X[0, :] = 3.0
    econ = np.zeros(n_rows - 1)
    Xo, _, _ = old_guard(X, econ)
    Xn, _, _ = new_guard(X, econ)
    assert np.array_equal(Xo[0], X[0]) and np.array_equal(Xn[0], X[0]), \
        "penalty row was modified"
    print("  penalty row untouched: OK")

    # 8. fuzz
    for trial in range(300):
        nr = int(rng.integers(2, 30))
        no = int(rng.integers(1, 9))
        X = rng.normal(size=(nr, no))
        X[rng.random(X.shape) < 0.4] = 0.0
        econ = np.abs(rng.normal(size=nr - 1))
        econ[rng.random(nr - 1) < 0.3] = 0.0
        check_silent(f"fuzz {trial}", X, econ)
    print("  300 random cases: OK")

    print("test_zero_error_guard OK")


def check_silent(name, X, econ_body):
    Xo, fired_o, rep_o = old_guard(X, econ_body)
    Xn, fired_n, rep_n = new_guard(X, econ_body)
    assert fired_o == fired_n, f"{name}: guard fired {fired_o} vs {fired_n}"
    assert np.array_equal(Xo, Xn), f"{name}: zeroed matrices differ"
    if rep_o is not None:
        for i in range(3):
            assert np.array_equal(np.sort(np.asarray(rep_o[i])),
                                  np.sort(np.asarray(rep_n[i]))), name


if __name__ == "__main__":
    main()
