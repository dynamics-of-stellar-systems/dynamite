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
    # bind the real NNLS methods so internal self.* calls resolve
    import types as _t

    for name in (
        "construct_nnls_matrix_and_rhs",
        "construct_adelie_matrix_and_rhs",
        "_prepare_kinematic_block",
        "apply_CR_cut",
        "_build_augmented_X",
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


if __name__ == "__main__":
    test_bitwise(np.float64)
    print("float64 OK")
    test_bitwise(np.float32)
    print("float32 OK")
    print("test_fused_assembly OK")
