"""Regression: bulk parsers must stay file-layout-synchronized under
streaming (subset-of-sets) reads.

Contract being pinned: ``ap_global`` must ALWAYS enumerate every aperture
stored in the file, in file order - ``n_pairs = norb * len(ap_global)`` is
what the record walk assumes. Streaming reads (``read_vel_histograms(
kin_sets=[i])``) therefore pass the FULL per-file list plus a boolean
``keep`` mask over it; only kept apertures scatter into ``velhist0``.

Root cause of the 2026-08-22 fused-stream chi2 corruption (8.1e9 instead of
2.9e6): read_orbit_base passed a SET-subset as ap_global, so the walk ran
over norb * |subset| pairs while the file held norb * |all 1d apertures|
records - every orbit after the first was read from shifted bytes, silently,
with no exception. Synthetic tests missed it because their stub replaced
read_orbit_base wholesale.

Run from the repo root:
    PYTHONPATH=. python dev_tests/test_streamed_parser_subset.py
"""

import os
import tempfile

import numpy as np

import test_hist_bulk_read as hbr
from dynamite.orblib import LegacyOrbitLibrary

# interleaved like omega Cen: 1d, 2d, 1d, 2d - so each FILE holds the
# apertures of TWO different kinematic sets, non-contiguously
KIN_SPECS = [(1, 2, 5), (2, 2, (3, 3)), (1, 3, 7), (2, 1, (5, 5))]
NORB = 4


def _alloc(kin_specs, norb):
    out = []
    for dim, n_ap, bins in kin_specs:
        if dim == 1:
            out.append(np.zeros((norb, bins, n_ap)))
        else:
            out.append(np.zeros((norb, bins[0], bins[1], n_ap)))
    return out


def _parse(lib, which_file, keep_set):
    """Parse one file of KIN_SPECS, keeping only kinematic set `keep_set`.

    Returns (velhist0, keep_mask_used)."""
    rng = np.random.default_rng(11)
    with tempfile.TemporaryDirectory() as td:
        losvd = os.path.join(td, "losvd.dat")
        pm = os.path.join(td, "pm.dat")
        truth, kin_idx_per_ap, idx_ap_reset, hist_bins, dims = hbr.build(
            KIN_SPECS, NORB, rng, losvd, pm
        )
        lib = lib or hbr._lib()
        ap_all = np.arange(len(kin_idx_per_ap))
        dim_per_ap = np.asarray(dims)[kin_idx_per_ap]
        fdim = 1 if which_file == "losvd" else 2
        ap_file = ap_all[dim_per_ap == fdim]
        keep = np.isin(kin_idx_per_ap[ap_file], [keep_set])
        velhist0 = [
            arr if i == keep_set else None
            for i, arr in enumerate(_alloc(KIN_SPECS, NORB))
        ]
        fname = losvd if fdim == 1 else pm
        if fdim == 1:
            lib._read_losvd_hist_vectorised(
                fname,
                NORB,
                ap_file,
                kin_idx_per_ap,
                idx_ap_reset,
                hist_bins,
                velhist0,
                keep=keep,
            )
        else:
            lib._read_pm_hist_vectorised(
                fname,
                NORB,
                ap_file,
                kin_idx_per_ap,
                idx_ap_reset,
                hist_bins,
                velhist0,
                keep=keep,
            )
        return velhist0[keep_set], truth[keep_set]


def test_subset_parse_matches_full():
    lib = LegacyOrbitLibrary.__new__(LegacyOrbitLibrary)
    for which_file, keep_set in (("losvd", 0), ("losvd", 2), ("pm", 1), ("pm", 3)):
        got, want = _parse(lib, which_file, keep_set)
        assert np.array_equal(got, want), (
            f"{which_file} set {keep_set}: streamed-subset parse differs "
            "from ground truth"
        )


if __name__ == "__main__":
    test_subset_parse_matches_full()
    print("test_streamed_parser_subset OK")
