"""The bulk histogram parsers must reproduce the per-record loop exactly.

These parsers read a binary file by walking Fortran record markers. A mistake
does not raise - it silently puts values in the wrong place, or reads the
wrong bytes - so this compares against a transcription of the loop from
``read_orbit_base`` over a range of library layouts, not just one.

Covered here:
  - 1d only, 2d only, and mixed libraries with SEVERAL sets of each, since a
    mixed library sends each parser a subset of apertures and the destination
    index depends on which kinematic set an aperture belongs to
  - sets interleaved as 1d, 2d, 1d, 2d (the omega Cen production layout), so
    the aperture subsets are not contiguous blocks
  - empty sentinel pairs, all-empty files, single orbit, single aperture
  - batch sizes from 1 upward, to split the scatter into many batches
  - a truncated file, which must raise rather than read past the end

Run this file.
"""
import logging
import os
import tempfile

import numpy as np
from scipy.io import FortranFile

from dynamite.orblib import LegacyOrbitLibrary

EMPTY_1D = np.array([8, -9], dtype=np.int32)
EMPTY_2D = np.array([8, 8, -9, -9], dtype=np.int32)


def _lib():
    """A LegacyOrbitLibrary with just enough state for the parsers."""
    lib = LegacyOrbitLibrary.__new__(LegacyOrbitLibrary)
    lib.mod_dir = 'test/'
    lib.logger = logging.getLogger('test_hist_bulk_read')
    lib.logger.addHandler(logging.NullHandler())
    return lib


def build(kin_specs, norb, rng, losvd_path, pm_path, empty_frac=0.5):
    """Write both files for a library and return the expected histograms.

    kin_specs : list of (dim, n_apertures, bins) - bins is an int for dim 1,
    a pair for dim 2, mirroring ``hist_bins`` in read_orbit_base.
    """
    dims = [d for d, _, _ in kin_specs]
    n_aps = [n for _, n, _ in kin_specs]
    hist_bins = [b for _, _, b in kin_specs]

    kin_idx_per_ap = np.concatenate(
        [np.zeros(n, dtype=int) + i for i, n in enumerate(n_aps)])
    idx_ap_reset = np.concatenate(([0], np.cumsum(n_aps)[:-1]))

    truth = []
    for kin, (dim, n_ap, bins) in enumerate(kin_specs):
        if dim == 1:
            truth.append(np.zeros((norb, bins, n_ap)))
        else:
            truth.append(np.zeros((norb, bins[0], bins[1], n_ap)))

    f1 = FortranFile(losvd_path, 'w')
    # read_orbit_base consumes one header record from the losvd file, and the
    # vectorised parser skips it too. The pm file has NO header.
    f1.write_record(np.array([0, 0], dtype=np.int32), np.array([0.0]))
    f2 = FortranFile(pm_path, 'w')
    try:
        for j in range(norb):
            for a_global, kin in enumerate(kin_idx_per_ap):
                a_local = a_global - idx_ap_reset[kin]
                bins = hist_bins[kin]
                if dims[kin] == 1:
                    c = (bins - 1) // 2
                    if rng.random() < empty_frac:
                        f1.write_record(EMPTY_1D)
                        continue
                    lo = int(rng.integers(-c, 1))
                    hi = int(rng.integers(lo, c + 1))
                    f1.write_record(np.array([lo, hi], dtype=np.int32))
                    vals = rng.random(hi - lo + 1)
                    f1.write_record(vals)
                    truth[kin][j, lo + c:hi + c + 1, a_local] = vals
                else:
                    c0, c1 = (bins[0] - 1) // 2, (bins[1] - 1) // 2
                    if rng.random() < empty_frac:
                        f2.write_record(EMPTY_2D)
                        continue
                    i0 = int(rng.integers(-c0, 1))
                    i1 = int(rng.integers(-c1, 1))
                    x0 = int(rng.integers(i0, c0 + 1))
                    x1 = int(rng.integers(i1, c1 + 1))
                    f2.write_record(np.array([i0, i1, x0, x1], dtype=np.int32))
                    block = rng.random((x0 - i0 + 1, x1 - i1 + 1))
                    f2.write_record(block.ravel(order='F'))  # Fortran order
                    truth[kin][j, i0 + c0:x0 + c0 + 1,
                               i1 + c1:x1 + c1 + 1, a_local] = block
    finally:
        f1.close()
        f2.close()
    return truth, kin_idx_per_ap, idx_ap_reset, hist_bins, dims


def reference_loop(losvd_path, pm_path, norb, kin_idx_per_ap, idx_ap_reset,
                   hist_bins, dims, out):
    """read_orbit_base's per-record loop, transcribed."""
    f1 = FortranFile(losvd_path, 'r')
    f1.read_record(np.int32, np.int32, float)          # header
    f2 = FortranFile(pm_path, 'r')
    try:
        for j in range(norb):
            for a_global, kin in enumerate(kin_idx_per_ap):
                a_local = a_global - idx_ap_reset[kin]
                bins = hist_bins[kin]
                if dims[kin] == 1:
                    lo, hi = f1.read_ints(np.int32)
                    if lo <= hi:
                        c = (bins - 1) // 2
                        out[kin][j, lo + c:hi + c + 1, a_local] = \
                            f1.read_reals(float)
                else:
                    i0, i1, x0, x1 = f2.read_ints(np.int32)
                    if i0 <= x0 and i1 <= x1:
                        c0, c1 = (bins[0] - 1) // 2, (bins[1] - 1) // 2
                        tmp = f2.read_reals(float).reshape(
                            (x0 - i0 + 1, x1 - i1 + 1), order='F')
                        out[kin][j, i0 + c0:x0 + c0 + 1,
                                 i1 + c1:x1 + c1 + 1, a_local] = tmp
    finally:
        f1.close()
        f2.close()


def bulk(losvd_path, pm_path, norb, kin_idx_per_ap, idx_ap_reset, hist_bins,
         dims, out, chunk):
    """The aperture-subset logic from read_orbit_base, plus both parsers."""
    lib = _lib()
    ap_all = np.arange(len(kin_idx_per_ap))
    dim_per_ap = np.asarray(dims)[kin_idx_per_ap]
    ap_1d, ap_2d = ap_all[dim_per_ap == 1], ap_all[dim_per_ap == 2]
    if ap_1d.size:
        lib._read_losvd_hist_vectorised(losvd_path, norb, ap_1d,
                                        kin_idx_per_ap, idx_ap_reset,
                                        hist_bins, out, chunk=chunk)
    if ap_2d.size:
        lib._read_pm_hist_vectorised(pm_path, norb, ap_2d, kin_idx_per_ap,
                                     idx_ap_reset, hist_bins, out, chunk=chunk)


LAYOUTS = {
    'omega Cen-like, interleaved 1d/2d':
        [(1, 3, 25), (2, 4, (5, 7)), (1, 2, 9), (2, 3, (3, 5))],
    '1d only, several sets':      [(1, 3, 7), (1, 4, 11), (1, 1, 5)],
    '2d only, several sets':      [(2, 2, (3, 5)), (2, 3, (7, 3))],
    '2d first, then 1d':          [(2, 2, (5, 5)), (1, 3, 9)],
    'single aperture each':       [(1, 1, 5), (2, 1, (3, 3))],
    'wide 2d, narrow 1d':         [(2, 5, (9, 11)), (1, 2, 3)],
}


def check(name, kin_specs, norb, empty_frac, chunks, rng):
    d = tempfile.mkdtemp()
    losvd, pm = os.path.join(d, 'l.dat'), os.path.join(d, 'p.dat')
    try:
        truth, kipa, reset, bins, dims = build(
            kin_specs, norb, rng, losvd, pm, empty_frac)

        ref = [np.zeros_like(t) for t in truth]
        reference_loop(losvd, pm, norb, kipa, reset, bins, dims, ref)
        for r, t in zip(ref, truth):
            assert np.array_equal(r, t), f'{name}: reference loop disagrees'

        for chunk in chunks:
            got = [np.zeros_like(t) for t in truth]
            bulk(losvd, pm, norb, kipa, reset, bins, dims, got, chunk)
            for kin, (g, r) in enumerate(zip(got, ref)):
                assert np.array_equal(g, r), \
                    f'{name}: set {kin} differs at chunk={chunk}'
        print(f'  ok  {name} (norb={norb}, empty_frac={empty_frac})')
    finally:
        for f in (losvd, pm):
            if os.path.exists(f):
                os.remove(f)
        os.rmdir(d)


def check_truncated():
    """A short file must raise ValueError, not read past the end."""
    rng = np.random.default_rng(7)
    d = tempfile.mkdtemp()
    losvd, pm = os.path.join(d, 'l.dat'), os.path.join(d, 'p.dat')
    try:
        specs = [(2, 3, (5, 5))]
        truth, kipa, reset, bins, dims = build(specs, 6, rng, losvd, pm, 0.3)
        with open(pm, 'r+b') as fh:      # lop off the last third
            fh.truncate(int(os.path.getsize(pm) * 0.66))
        out = [np.zeros_like(t) for t in truth]
        try:
            bulk(losvd, pm, 6, kipa, reset, bins, dims, out, 1000)
        except ValueError as e:
            assert 'truncated' in str(e), f'unexpected message: {e}'
            print('  ok  truncated pm file raises ValueError')
        else:
            raise AssertionError('truncated file did not raise')
    finally:
        for f in (losvd, pm):
            if os.path.exists(f):
                os.remove(f)
        os.rmdir(d)


def demo():
    rng = np.random.default_rng(0)
    print('layouts:')
    for name, specs in LAYOUTS.items():
        for norb, empty in ((5, 0.5), (1, 0.0), (3, 1.0)):
            check(name, specs, norb, empty, (1, 7, 10**9), rng)
    print('edge cases:')
    check_truncated()
    print('\nbulk histogram read == per-record loop, OK')


if __name__ == '__main__':
    demo()
