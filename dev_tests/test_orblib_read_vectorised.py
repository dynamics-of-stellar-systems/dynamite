#!/usr/bin/env python3
"""Check the vectorised LOSVD read against the original per-record loop.

``LegacyOrbitLibrary.read_orbit_base`` parses ``{fileroot}_losvd_hist.dat`` in
bulk (see ``_read_losvd_hist_vectorised``) instead of doing two
``scipy.io.FortranFile`` record reads per (orbit, aperture) pair. The two must
agree exactly, so this reads a model's orbit library twice -- once with the
vectorised parser monkeypatched back to the original loop, once with the real
thing -- and compares the LOSVDs, intrinsic masses and projected masses.

Needs a model tree with an already-integrated orbit library; it does not
integrate anything.

Usage:
    python test_orblib_read_vectorised.py <config.yaml> [--dir <workdir>]
"""

import argparse
import os
import sys
import time

import numpy as np
from scipy.io import FortranFile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import dynamite as dyn                                    # noqa: E402
from dynamite.orblib import LegacyOrbitLibrary            # noqa: E402


def reference_loop(self, fname, norb, kin_idx_per_ap, idx_ap_reset,
                   hist_bins, velhist0, chunk=None):
    """The original read loop, kept as the reference implementation."""
    fort_file = FortranFile(fname, 'r')
    _ = fort_file.read_record(np.int32, np.int32, float)   # header record
    for j in range(norb):
        for i_ap, kin_idx in enumerate(kin_idx_per_ap):
            i_ap0 = i_ap - idx_ap_reset[kin_idx]
            ivmin, ivmax = fort_file.read_ints(np.int32)
            if ivmin <= ivmax:
                nv0 = (hist_bins[kin_idx] - 1) // 2
                velhist0[kin_idx][j, ivmin + nv0:ivmax + nv0 + 1, i_ap0] = \
                    fort_file.read_reals(float)
    fort_file.close()


def read_once(config_file, row, use_reference):
    """Read a model's orbit library; return its arrays and the elapsed time."""
    real = LegacyOrbitLibrary._read_losvd_hist_vectorised
    if use_reference:
        LegacyOrbitLibrary._read_losvd_hist_vectorised = reference_loop
    try:
        c = dyn.config_reader.Configuration(config_file, reset_logging=True,
                                            reset_existing_output=False)
        if len(c.all_models.table) > row:
            mod = c.all_models.get_model_from_row(row)
        else:
            # a tree with an integrated orblib but never run through
            # ModelIterator has an empty all_models.ecsv, so build the parset
            # straight from the (fixed) parameter space instead
            from astropy import table as aptable
            vals = {p.name: [p.par_value] for p in c.parspace}
            parset = aptable.Table(
                {n: vals[n] for n in c.parspace.par_names})[0]
            mod = dyn.model.Model(config=c, parset=parset)
        orblib = mod.get_orblib()
        t_0 = time.perf_counter()
        orblib.read_vel_histograms()
        elapsed = time.perf_counter() - t_0
        arrays = dict(
            losvd=[np.asarray(h.y).copy() for h in orblib.vel_histograms],
            intrinsic=np.asarray(orblib.intrinsic_masses).copy(),
            projected=[np.asarray(p).copy()
                       for p in orblib.projected_masses])
        return arrays, elapsed
    finally:
        LegacyOrbitLibrary._read_losvd_hist_vectorised = real


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    parser.add_argument('--dir', default=None,
                        help='change to this directory first')
    parser.add_argument('--row', type=int, default=0,
                        help='all_models row index')
    args = parser.parse_args()
    if args.dir:
        os.chdir(args.dir)

    print(f'dynamite {dyn.__version__} from {dyn.__path__[0]}', flush=True)

    print('[1/2] original per-record loop ...', flush=True)
    ref, t_ref = read_once(args.config, args.row, use_reference=True)
    print(f'      {t_ref:.1f}s', flush=True)

    print('[2/2] vectorised parser ...', flush=True)
    new, t_new = read_once(args.config, args.row, use_reference=False)
    print(f'      {t_new:.1f}s', flush=True)

    ok = True
    for i, (a, b) in enumerate(zip(ref['losvd'], new['losvd'])):
        same = a.shape == b.shape and np.array_equal(a, b)
        ok &= same
        print(f'  losvd set {i}: shape {a.shape} '
              f'nonzero {np.count_nonzero(a)} identical={same}')
    for key in ('intrinsic', 'projected'):
        a, b = ref[key], new[key]
        if isinstance(a, list):
            same = (len(a) == len(b)
                    and all(np.array_equal(x, y) for x, y in zip(a, b)))
        else:
            same = np.array_equal(a, b)
        ok &= same
        print(f'  {key:10s}: identical={same}')

    print(f'\nbit-identical : {ok}')
    print(f'read time     : {t_ref:.1f}s -> {t_new:.1f}s '
          f'({t_ref / max(t_new, 1e-9):.1f}x)')
    if not ok:
        print('FAILED: the vectorised parser does not reproduce the loop.')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
