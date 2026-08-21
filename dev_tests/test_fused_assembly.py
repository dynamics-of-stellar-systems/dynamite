"""Bitwise equivalence of the fused adelie constructor vs the two-step one.

    construct_adelie_matrix_and_rhs(orblib)
        == _build_augmented_X(*construct_nnls_matrix_and_rhs(orblib))

Drives the REAL methods on a minimal fake system / kinematic sets / orblib
(real Histogram objects; transforms return non-contiguous moveaxis'd views so
the reshape-view asserts are genuinely exercised).

Note on ordering: the reference call runs first and its in-place edge-bin
zeroing (1D sets) persists on the shared histograms; the fused call then sees
pre-zeroed histograms. Zeroing is idempotent, so both paths still see
identical values and the comparison stays exact.

Run from the repo root:
    PYTHONPATH=. python dev_tests/test_fused_assembly.py
"""

import logging
import types

import numpy as np

from dynamite import kinematics as dyn_kin
from dynamite.weight_solvers import NNLS

N_ORBS = 7


class FakeKins1D:
    """GH-like: 3 spatial bins, 4 constraints per orbit per aperture."""

    kind = "gh"
    n_spatial_bins = 3
    constraint_shape = (4,)

    def get_observed_values_and_uncertainties(self, settings):
        obs = np.arange(self.n_spatial_bins * 4, dtype=float).reshape(3, 4)
        err = np.abs(obs) * 0.01 + 0.01
        return obs, err

    def transform_orblib_to_observables(self, hist, settings):
        # non-contiguous view carrying ALL apertures per orbit, like the real
        # implementations: (n_orbs, aps, c) from (n_orbs, c', aps)
        return np.moveaxis(hist.y[:, :4, :], 1, 2)


class FakeKins2D:
    """PM-like: 2 spatial bins, 16 constraints per aperture."""

    kind = "pm"
    n_spatial_bins = 2
    constraint_shape = (16,)

    def get_observed_values_and_uncertainties(self, settings):
        obs = np.arange(self.n_spatial_bins * 16, dtype=float).reshape(2, 16)
        err = np.abs(obs) * 0.02 + 0.02
        return obs, err

    def transform_orblib_to_observables(self, hist, settings):
        # moveaxis'd non-contiguous view: (n_orbs, aps, 4, 4)
        return np.moveaxis(hist.y[:, :4, :4, :], 3, 1)


def _make_fake_nnls(dtype):
    mge_intrinsic = np.array([3e5, 2e5, 1e5, 0.5e5])
    mge_projected = np.array([4e5, 3e5, 2e5, 0.7e5, 0.3e5])
    fake = types.SimpleNamespace()
    fake.nnls_dtype = dtype
    fake.adelie_mu = 1e7
    fake.adelie_tol = 1e-10
    fake.adelie_gap_tol = 1e-10
    fake.adelie_alm_iters = 5
    fake.CRcut = False
    fake.settings = {}
    fake.n_intrinsic = 4
    fake.n_apertures = 5
    fake.n_mass_constraints = 1 + fake.n_intrinsic + fake.n_apertures
    fake.intrinsic_masses = mge_intrinsic
    fake.intrinsic_mass_error = 0.05
    fake.projected_masses = mge_projected
    fake.projected_mass_error = 0.1
    fake.total_mass = 1e6
    fake.total_mass_error = max(abs(1.0 - 1e6), 1e-8)
    fake.direc_with_ml = "fake"
    fake.logger = logging.getLogger("test_fused_assembly")
    fake.system = types.SimpleNamespace(
        get_unique_triaxial_visible_component=lambda: types.SimpleNamespace(
            kinematic_data=[FakeKins1D(), FakeKins2D()]
        )
    )
    # bind the real NNLS methods so internal self.* calls resolve.
    # Staticmethods must be attached as plain functions - wrapping them in
    # MethodType would prepend the fake instance as an extra first argument.
    import types as _t

    for name in (
        "_build_augmented_X",
        "kkt_violation",
        "kkt_violation_augmented",
    ):
        setattr(fake, name, getattr(NNLS, name))
    for name in (
        "construct_nnls_matrix_and_rhs",
        "construct_adelie_matrix_and_rhs",
        "_prepare_kinematic_block",
        "apply_CR_cut",
        "solve_adelie_alm",
    ):
        setattr(fake, name, _t.MethodType(getattr(NNLS, name), fake))
    return fake


def _make_fake_orblib():
    rng = np.random.default_rng(5)
    # 1D hist: y (n_orbs, nv=6, aps=3); edges get zeroed destructively by the
    # constructors, so start them nonzero to make that zeroing observable
    y1 = rng.random((N_ORBS, 6, 3))
    hist1 = dyn_kin.Histogram(xedg=np.linspace(-3.0, 3.0, 7), y=y1)
    # 2D hist: y (n_orbs, nv0=4, nv1=5, aps=2)
    y2 = rng.random((N_ORBS, 4, 5, 2))
    hist2 = dyn_kin.Histogram2D(
        xedg=(np.linspace(-2.0, 2.0, 5), np.linspace(-2.5, 2.5, 6)), y=y2
    )
    orblib = types.SimpleNamespace()
    orblib.n_orbs = N_ORBS
    orblib.vel_histograms = [hist1, hist2]
    orblib.intrinsic_masses = rng.random((N_ORBS, 4))
    orblib.projected_masses = [
        rng.random((N_ORBS, 3)) * 1e5,
        rng.random((N_ORBS, 2)) * 1e5,
    ]
    return orblib


def test_bitwise(dtype):
    mu = 1e7
    sqrt_mu = np.sqrt(mu)
    fake = _make_fake_nnls(dtype)
    orblib = _make_fake_orblib()

    A_ref, b_ref = NNLS.construct_nnls_matrix_and_rhs(fake, orblib)
    X_ref, cn_ref, y_ref = NNLS._build_augmented_X(A_ref[1:], b_ref[1:], sqrt_mu, dtype)

    prob = NNLS.construct_adelie_matrix_and_rhs(fake, orblib)

    assert prob.X.shape == X_ref.shape
    assert np.array_equal(prob.X, X_ref), f"X differs ({dtype})"
    assert np.array_equal(prob.col_norm, cn_ref), f"col_norm differs ({dtype})"
    assert np.array_equal(prob.y, y_ref), f"y differs ({dtype})"
    assert np.array_equal(prob.row0_vec, A_ref[0]), f"row0_vec != A[0] ({dtype})"
    assert prob.b0 == float(b_ref[0])


def test_alm_end_to_end(dtype):
    """The full new pipeline on synthetic data: fused problem -> ALM ->
    surrogate chi2_vector / KKT, compared against the A-based forms."""
    from dynamite.weight_solvers import chi2_vector_from_residuals

    fake = _make_fake_nnls(dtype)
    orblib = _make_fake_orblib()
    A_ref, b_ref = NNLS.construct_nnls_matrix_and_rhs(fake, orblib)
    prob = NNLS.construct_adelie_matrix_and_rhs(fake, orblib)

    w, resid_full = NNLS.solve_adelie_alm(fake, prob)
    assert w.shape == (prob.X.shape[1],)
    assert np.all(w >= 0)
    assert resid_full.shape == A_ref.shape[0:1]

    row0_resid = float(prob.row0_vec @ w) - prob.b0
    chi2_new = float(np.sum(chi2_vector_from_residuals(resid_full, row0_resid**2)))
    chi2_ref = float(np.sum((A_ref @ w - b_ref) ** 2))
    # float64: round-off only. float32: accumulated gemv/accumulation
    # rounding of the scaled matrix shows up at ~1e-11 relative.
    chi2_tol = 1e-11 if dtype == np.float64 else 1e-7
    assert np.isclose(chi2_new, chi2_ref, rtol=chi2_tol), (chi2_new, chi2_ref)

    kkt_new = NNLS.kkt_violation_augmented(
        prob.row0_vec,
        prob.b0,
        prob.X,
        prob.col_norm,
        resid_full,
        w,
        fake.adelie_mu,
    )
    kkt_ref = NNLS.kkt_violation(A_ref.astype(np.float64), b_ref.astype(np.float64), w)
    assert 0.0 <= kkt_new[0] <= 1.0
    # The scaled value is a ratio of cancelling terms, so float32 agrees far
    # more loosely than chi2; raw agrees tightly in both dtypes.
    kkt_tol = 1e-8 if dtype == np.float64 else 1e-3
    assert np.isclose(kkt_new[0], kkt_ref[0], rtol=kkt_tol), (kkt_new, kkt_ref)
    assert np.isclose(kkt_new[1], kkt_ref[1], rtol=1e-6), (kkt_new, kkt_ref)


if __name__ == "__main__":
    test_bitwise(np.float64)
    print("float64 OK")
    test_bitwise(np.float32)
    print("float32 OK")
    test_alm_end_to_end(np.float64)
    test_alm_end_to_end(np.float32)
    print("alm end-to-end OK")
    print("test_fused_assembly OK")
