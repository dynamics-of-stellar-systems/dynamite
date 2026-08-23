"""read_orbit_base's per-file gates: need_1d / need_2d and the keep masks.

Streaming reads one kinematic set at a time, and 1d and 2d histograms live in
separate files, so a call that wants only 1d sets must not decompress and walk
the pm file. Two things have to hold, and neither is covered by the parser
tests (which call the parsers directly, with read_orbit_base stubbed out):

1. With kin_sets=None (requested_sets = every set) the gates must reduce
   EXACTLY to the old `any(i == N for i in hist_dim)` conditions, so the
   non-streaming path is unchanged.
2. need_2d False must imply keep_2d.any() False (and likewise for 1d).
   `tmpfname_pm` only exists when need_2d opened it, but it is referenced
   under the keep gate - if the two could disagree, that is an unbound
   local at production scale only.

Run: python test_read_gate_logic.py
"""

import itertools

import numpy as np


def gates(hist_dim, requested_sets, n_apertures):
    """The decisions read_orbit_base makes, transcribed."""
    need_1d = any(hist_dim[si] == 1 for si in requested_sets)
    need_2d = any(hist_dim[si] == 2 for si in requested_sets)

    kin_idx_per_ap = np.array(
        [si for si, n in enumerate(n_apertures) for _ in range(n)]
    )
    ap_all = np.arange(len(kin_idx_per_ap))
    hist_dim_per_ap = np.asarray(hist_dim)[kin_idx_per_ap]
    ap_1d = ap_all[hist_dim_per_ap == 1]
    ap_2d = ap_all[hist_dim_per_ap == 2]
    keep_1d = np.isin(kin_idx_per_ap[ap_1d], requested_sets)
    keep_2d = np.isin(kin_idx_per_ap[ap_2d], requested_sets)
    return need_1d, need_2d, ap_1d, ap_2d, keep_1d, keep_2d


def main():
    rng = np.random.default_rng(20260823)
    n_checked = 0

    # every layout of up to 4 kinematic sets, every subset of them
    for n_kins in range(1, 5):
        for hist_dim in itertools.product([1, 2], repeat=n_kins):
            hist_dim = list(hist_dim)
            n_apertures = [int(rng.integers(1, 5)) for _ in range(n_kins)]

            # 1. kin_sets=None must reproduce the old conditions exactly
            all_sets = list(range(n_kins))
            need_1d, need_2d, *_ = gates(hist_dim, all_sets, n_apertures)
            assert need_1d == any(i == 1 for i in hist_dim), (hist_dim, "1d")
            assert need_2d == any(i == 2 for i in hist_dim), (hist_dim, "2d")

            # 2. over every non-empty subset of sets
            for r in range(1, n_kins + 1):
                for subset in itertools.combinations(all_sets, r):
                    requested = sorted(subset)
                    (need_1d, need_2d, ap_1d, ap_2d,
                     keep_1d, keep_2d) = gates(hist_dim, requested, n_apertures)

                    # the load-bearing coupling: the file is opened iff the
                    # keep mask can select anything from it
                    assert need_1d == bool(ap_1d.size and keep_1d.any()), \
                        (hist_dim, requested, "1d gate/keep disagree")
                    assert need_2d == bool(ap_2d.size and keep_2d.any()), \
                        (hist_dim, requested, "2d gate/keep disagree")

                    # and the kept apertures are exactly the requested sets'
                    kin_idx_per_ap = np.array(
                        [si for si, n in enumerate(n_apertures)
                         for _ in range(n)]
                    )
                    for aps, keep, dim in ((ap_1d, keep_1d, 1),
                                           (ap_2d, keep_2d, 2)):
                        want = {a for a in aps
                                if kin_idx_per_ap[a] in requested}
                        got = set(np.asarray(aps)[keep])
                        assert want == got, (hist_dim, requested, dim)
                    n_checked += 1

    print(f"  {n_checked} (layout, subset) combinations: OK")
    print("test_read_gate_logic OK")


if __name__ == "__main__":
    main()
