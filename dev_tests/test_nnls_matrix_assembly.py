"""Check the pre-sized, F-ordered assembly reproduces the old vstack version.

The real construct_nnls_matrix_and_rhs needs a stored orbit library (see
_real_orblib_check.py for the end-to-end version); this mirrors just its
assembly arithmetic on synthetic blocks. Run: python this file.
"""
import numpy as np


def _blocks(rng, n_orbs, n_mass, block_shapes):
    mass = rng.normal(size=(n_mass, n_orbs))
    con0 = rng.normal(size=n_mass)
    econ0 = np.abs(rng.normal(size=n_mass)) + 0.1
    kin = []
    for shape in block_shapes:
        # as transform_orblib_to_observables returns it: a moveaxis'd,
        # non-contiguous view of shape (n_orbs, ...)
        raw = rng.normal(size=(n_orbs,) + shape)
        orb_kins = np.moveaxis(raw, -1, 1) if len(shape) > 1 else raw
        n_c = orb_kins.size // n_orbs
        kin.append((orb_kins, rng.normal(size=n_c),
                    np.abs(rng.normal(size=n_c)) + 0.1))
    return mass, con0, econ0, kin


def old_way(mass, con0, econ0, kin, dtype):
    con, econ = con0.astype(dtype), econ0.astype(dtype)
    orbmat = mass.astype(dtype)
    n_orbs = mass.shape[1]
    for orb_kins, obs_kins, obs_err in kin:
        con = np.concatenate((con, obs_kins.astype(dtype)))
        econ = np.concatenate((econ, obs_err.astype(dtype)))
        flat = np.reshape(orb_kins, (n_orbs, -1)).astype(dtype)
        orbmat = np.vstack((orbmat, flat.T))
    return con, econ, orbmat


def new_way(mass, con0, econ0, kin, dtype):
    n_mass, n_orbs = mass.shape
    n_rows = n_mass + sum(o.size for _, o, _ in kin)
    con = np.zeros(n_rows, dtype=dtype)
    econ = np.zeros(n_rows, dtype=dtype)
    orbmat = np.zeros((n_rows, n_orbs), dtype=dtype, order='F')
    con[:n_mass], econ[:n_mass], orbmat[:n_mass] = con0, econ0, mass
    i = n_mass
    for orb_kins, obs_kins, obs_err in kin:
        j = i + obs_kins.size
        assert orb_kins.size // n_orbs == obs_kins.size
        con[i:j], econ[i:j] = obs_kins, obs_err
        dest = orbmat[i:j, :].T.reshape(orb_kins.shape)
        assert np.shares_memory(dest, orbmat), 'block write is not a view'
        dest[...] = orb_kins
        i = j
    return con, econ, orbmat


def demo():
    rng = np.random.default_rng(0)
    shapes = ((11,), (4, 3), (2, 5, 3))   # 1D GH-like and 2D/3D PM-like blocks
    for dtype in (np.float64, np.float32):
        args = _blocks(rng, n_orbs=37, n_mass=5, block_shapes=shapes)
        c1, e1, m1 = old_way(*args, dtype)
        c2, e2, m2 = new_way(*args, dtype)
        assert m1.dtype == m2.dtype == dtype
        assert np.array_equal(c1, c2) and np.array_equal(e1, e2)
        assert np.array_equal(m1, m2), 'destination-view write disagrees'
        rhs1 = np.zeros_like(c1)
        np.divide(c1, e1, out=rhs1, where=e1 != 0)
        t1, t2 = m1.T.copy(), m2.T.copy()
        np.divide(t1, e1, out=t1, where=e1 != 0)
        np.divide(t2, e2, out=t2, where=e2 != 0)
        assert np.array_equal(t1, t2)

    print("nnls matrix assembly: F-ordered destination write == vstack, OK")


if __name__ == "__main__":
    demo()
