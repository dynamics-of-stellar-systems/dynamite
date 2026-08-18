"""The bulk 2d PM histogram parser must reproduce the per-record loop.

Writes a small pm_hist file with scipy's FortranFile - the same record framing
the Fortran code emits - then fills the destination array both ways and
compares. Run this file.
"""
import os
import tempfile

import numpy as np
from scipy.io import FortranFile

from dynamite.orblib import LegacyOrbitLibrary


def write_file(path, norb, n_ap, bins, rng, empty_frac=0.6):
    """Write records exactly as read_orbit_base's hist_dim == 2 branch expects."""
    truth = np.zeros((norb, bins[0], bins[1], n_ap))
    c0, c1 = (bins[0] - 1) // 2, (bins[1] - 1) // 2
    with FortranFile(path, 'w') as f:
        for j in range(norb):
            for a in range(n_ap):
                if rng.random() < empty_frac:            # empty sentinel
                    f.write_record(np.array([8, 8, -9, -9], dtype=np.int32))
                    continue
                i0 = int(rng.integers(-c0, 1))
                i1 = int(rng.integers(-c1, 1))
                x0 = int(rng.integers(i0, c0 + 1))
                x1 = int(rng.integers(i1, c1 + 1))
                f.write_record(np.array([i0, i1, x0, x1], dtype=np.int32))
                block = rng.random((x0 - i0 + 1, x1 - i1 + 1))
                # Fortran order on the wire, as the reader's reshape assumes
                f.write_record(block.ravel(order='F'))
                truth[j, i0 + c0:x0 + c0 + 1, i1 + c1:x1 + c1 + 1, a] = block
    return truth


def reference_loop(path, norb, n_ap, bins, out):
    """The original per-record read, transcribed."""
    c0, c1 = (bins[0] - 1) // 2, (bins[1] - 1) // 2
    with FortranFile(path, 'r') as f:
        for j in range(norb):
            for a in range(n_ap):
                i0, i1, x0, x1 = f.read_ints(np.int32)
                if i0 <= x0 and i1 <= x1:
                    tmp = f.read_reals(float)
                    tmp = tmp.reshape((x0 - i0 + 1, x1 - i1 + 1), order='F')
                    out[j, i0 + c0:x0 + c0 + 1, i1 + c1:x1 + c1 + 1, a] = tmp


def demo():
    rng = np.random.default_rng(0)
    for norb, n_ap, bins in ((7, 5, (5, 7)), (3, 1, (3, 3)), (4, 6, (9, 5))):
        fd, path = tempfile.mkstemp(suffix='_pm_hist.dat')
        os.close(fd)
        try:
            truth = write_file(path, norb, n_ap, bins, rng)

            ref = np.zeros_like(truth)
            reference_loop(path, norb, n_ap, bins, ref)
            assert np.array_equal(ref, truth), 'the reference loop is wrong'

            got = np.zeros_like(truth)
            lib = LegacyOrbitLibrary.__new__(LegacyOrbitLibrary)
            lib.mod_dir = 'test'
            lib._read_pm_hist_vectorised(
                path, norb,
                ap_global=np.arange(n_ap),
                kin_idx_per_ap=np.zeros(n_ap, dtype=int),
                idx_ap_reset=np.array([0]),
                hist_bins=[np.array(bins)],
                velhist0=[got],
                chunk=17)                       # tiny, to force many batches
            assert np.array_equal(got, ref), f'bulk read differs for {bins}'
        finally:
            os.remove(path)
    print("bulk 2d PM read == per-record loop, OK")


if __name__ == '__main__':
    demo()
