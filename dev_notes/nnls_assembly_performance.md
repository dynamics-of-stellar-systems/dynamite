# NNLS matrix assembly: memory and layout

Notes from profiling `NNLS.construct_nnls_matrix_and_rhs` and the adelie ALM
path against a real orbit library, 2026-08-18.

## Setup

`NGC5139_veldist_combined_bigger` (MUSE BayesLOSVD + HST proper motions,
3840 orbits) in `~/research/omegacen/dynamite_dataprep`. The assembled matrix
is 322403 x 3840 float64 = 9.9 GiB. Production omega Cen is ~45000 orbits and
~371212 constraints = ~125 GiB, i.e. about 12x this.

Reproduce with `dev_tests/_real_orblib_check.py <out.npz>`, run from the
dataprep directory with the `main` conda env. It profiles the assembly and
dumps `A`/`b`, so two git revisions can be compared with `np.array_equal`.

## Result

`construct_nnls_matrix_and_rhs`: **17.2s -> 3.4s**, matrix bit-identical.

Three changes, all layout/allocation, none numerical:

1. **Pre-size `con`/`econ`/`orbmat`** instead of growing `orbmat` with
   `np.vstack` per kinematic set. The row count per set is
   `get_observed_values_and_uncertainties(...)[0].size`, which depends only on
   the kinematic data, so a cheap pre-pass can size everything. Growing peaked
   at ~2.2x the matrix (old + new + block + its astype copy); pre-sizing peaks
   at ~1.4x. Time saved is minor (~1.5s); this one is about peak RSS.

2. **Allocate `orbmat` Fortran-ordered.** Blocks arrive as
   `(n_orbs, n_constraints)` and are written in transposed. Per 2.3 GiB block:

   | | C-ordered | F-ordered |
   |---|---|---|
   | block write | 1.00 s | 0.08 s |
   | divide by `econ` | 0.12 s | 0.08 s |
   | build adelie's F-contiguous `X` | 2.38 s | 0.33 s |
   | `A @ w` (chi2) | 0.04 s | 0.07 s |

3. **Reshape the destination, not the source.** `orb_kins` is a
   `moveaxis`/`swapaxes` view, so `np.reshape(orb_kins, (n_orbs, -1))` cannot
   be a view and copies the entire block (2.1s, and a full-size temporary -
   ~100 GiB at omega Cen scale). Reshaping the destination slice is a view, so
   the write needs no temporary: 1.33s/2314 MiB peak -> 1.03s/0 MiB peak.

   `ndarray.reshape` silently returns a **copy** when it cannot return a view,
   which as an assignment target is a silent no-op leaving the block at zero.
   The `np.shares_memory` assert in the code exists to make that loud. The
   `.shape =` setter raises natively but is deprecated in numpy 2.5, and the
   package supports numpy>=1.26.

## Also fixed, in the solver path

- `A_normalized`/`b_normalized` were built unconditionally but are unused on
  the adelie path, holding a second full copy of `A` for the whole run.
- `A.astype(np.float64)` in the ALM loop is a no-op cast at the default
  `nnls_dtype`, but `.astype` copies anyway - a full copy of `A` per iteration.
  Now `copy=False`.
- Building `X` went `np.vstack` -> `X / col_norm` -> `np.asfortranarray`, three
  full-size allocations; it now fills one F-ordered buffer in place. With `A`
  alive that peak was ~4x the matrix (~500 GiB at omega Cen scale), plus the
  unused `A_normalized` on top.

## Things checked and deliberately not changed

- **`Histogram2D.rebin_orblib_to_observations`** (the PM path) looks terrible -
  nested Python loops recomputing `overlapping_bins` inside the outer loop -
  and the overlap weights are exactly separable, so it collapses to two
  tensordots. Verified to 3e-16, but only 1.0-1.7x faster at 20x20 to 60x60
  velocity bins, because the overlaps stay local so the loop is effectively
  linear in bin count. Not worth the risk. In the test config it does not even
  run (`No rebinning necessary`). Revisit only for much finer PM grids.
- **`GaussHermite.get_gh_expansion_coefficients`**: `optimize=True` is only
  ~1.5x and the largest intermediate is ~108 MB with no orbit axis, so it is
  safe but not important. Applied anyway since it is one word.

## Where the time actually goes

Earlier guesses in this area (mine included) blamed `np.vstack` and then
`transform_orblib_to_observables` for a reported hour-long assembly phase.
Neither is right. On this library the transforms are noise - LOSVD rebin
0.13s, PM transform 0.002s - and the cost was in the plain assignment
statements, i.e. memory layout. Scaling 17s by 12x puts omega Cen assembly at
3-4 minutes even before these changes, so **whatever took an hour was not this
function**. `orblib.read_vel_histograms()` takes 77s here, 20x the assembly;
that is the next thing to profile.

## Not yet validated

The ALM loop itself has not been re-run end to end on real data, so the
predicted knock-on win from F-ordered `A` in `solve_adelie_alm` is inferred
from the microbenchmark above, not measured in situ.

## Still open

- No full `solve()` has been profiled on real data here. The solvers are slow
  enough on this laptop to look hung; they run fine on the compute cluster, so
  profile the solve there rather than chasing it locally.
- `nnls_solver: "cvxopt"` builds `P = A_normalized.T @ A_normalized`, an
  (n_orbs, n_orbs) Gram matrix. At omega Cen's 45000 orbits that is 16 GB and
  ~7.5e14 flops. Unusable at production scale; fine for small problems.
- The configs here set `VECLIB_MAXIMUM_THREADS=1` to dodge an Accelerate
  SIGSEGV in `scipy.optimize.nnls`. That also serialises every BLAS call in
  the solve, which is worth remembering when timing the scipy path.

## The ALM chi2 shortcut (2026-08-18, follow-up)

With adelie as the production solver for omega Cen, the per-iteration cost of
the ALM loop is what matters. Each multiplier update recomputed

    chi2 = sum((A @ w - b)**2)

which is a full pass over the matrix - 125 GiB read for omega Cen - repeated
up to `adelie_alm_iters` times, on top of the solve itself.

It is not needed. `X` is built as `A[1:]` with unit-L2 column scaling and
`w = beta / col_norm`, so `X[1:] @ beta == A[1:] @ w` exactly, and adelie's
returned state already carries `resid = y - X @ beta` (verified: equal to
`y - X @ beta` to 2e-15, for uniform and non-uniform observation weights
alike; `state.loss` is `0.5 * sum(weights * resid**2)`). Only row 0 of A needs
evaluating, because `X` replaces it with the ALM penalty row, and that is one
dot product over the orbits:

    chi2 = (A[0] @ w - b[0])**2 + resid[1:] @ resid[1:]

Per-iteration passes over A therefore go from one to zero. Outside the loop A
is still read three times per solve: the final `chi2_vector`, and the two
matvecs inside `kkt_violation`.

What that is worth, measured on the real 9.2 GiB matrix (F-ordered, as
`solve_adelie_alm` now receives it):

    old chi2 line, A @ w      2618 ms   (3.5 GiB/s effective on this laptop)
    new chi2 line, resid dot     0.04 ms

`adelie_alm_iters` defaults to 200 and the real subsampled solve used all 200,
so this is ~520s per model here. The cost is one full read of A per iteration,
so it scales with the matrix: omega Cen's ~125 GiB matrix means ~25 TiB of
reads per model, which is 8 to 40 minutes depending on the node's achievable
memory bandwidth (~50 GiB/s to ~10 GiB/s). The laptop's 3.5 GiB/s figure
should NOT be extrapolated directly - a production node is much wider.

Either way it is pure overhead outside the solver, and it is now gone.

Verified in `dev_tests/test_alm_chi2_from_resid.py`, which drives the real
`adelie.solver.bvls` through the same warm-started ALM loop and compares both
chi2 forms every iteration (agreement 0 to 3e-15 relative).

Also verified end to end by `dev_tests/_real_alm_chi2_check.py`, which runs
the real `solve_adelie_alm` on a row-subsampled real matrix and recomputes
chi2 directly from the returned weights. Two runs, both 200 ALM iterations:

    rows    reported (from resid)   direct A@w - b      final gap
    8000    16910.98                16910.9799656940    3.5e-09
    40000   91786.1095              91786.109536        7.3e-08

Exact to the four decimal places the solver logs in both cases. Not yet run at
full omega Cen scale, where the matrix has ~40x more rows than the larger of
these.

(An earlier version of this note claimed these runs died on the laptop. They
did not - they were killed by the way they were launched. They complete.)

## Possible follow-up: dropping A during the solve

After this change A is used inside `solve_adelie_alm` only for `A[0] @ w`
(a single row). If `chi2_vector` and `kkt_violation` were also expressed
through `X` and `col_norm`, A could be released once X is built, halving the
resident matrix count from two to one (~125 GiB). Not attempted: it changes
what `kkt_violation` is handed, and that function is the only optimality
certificate in the code.

## The next win: the bulk read is disabled by proper motions (2026-08-18)

With assembly and the ALM loop fixed, `orblib.read_vel_histograms()` is what
is left, and it is far slower than the data justifies. Profiled on the
NGC5139 combined library:

    read_vel_histograms: 163.2s, 119 million python calls

      ncalls    tottime  function
     5134940     71.05   numpy.fromfile
    10269878     24.38   _io.BufferedReader.read
     5134936     15.78   scipy/io/_fortran.py:170(read_record)
    10269872      9.09   scipy/io/_fortran.py:127(_read_size)
    10269872      7.81   numpy.frombuffer
           2     10.50   orblib.py:1542(read_orbit_base)
           2      4.95   orblib.py:1848(combine_and_mirror_orblibs)

That is 5.13 million individual scipy FortranFile record reads at ~32 us
each, roughly 111s of the 163s, against a few hundred MB of actual data
(the two qgrid files decompress to 157 MB each). The read is overhead-bound
by about two orders of magnitude, not I/O-bound.

The cause is a one-line gate in `read_orbit_base`:

    vectorised = (not legacy_file and not return_intrinsic_moments
                  and all(i == 1 for i in hist_dim))

`_read_losvd_hist_vectorised` - the bulk parser behind the earlier 24x LOSVD
read speedup - only runs when the library contains **nothing but 1D LOSVD
histograms**. Any proper-motion set makes `hist_dim` contain a 2, the gate
fails, and the whole library falls back to the per-(orbit, aperture) Python
loop. The code comment there already anticipates this ("tens of millions of
python-level calls").

So the optimisation is switched off precisely for the omega Cen production
configuration, which is MUSE LOSVDs plus two proper-motion sets. The 1405
apertures of the PM set dominate the call count.

The fix is to extend the bulk parser to 2D histograms, or failing that to
parse the 1D sets in bulk and loop only over the 2D ones. Not attempted here.

## Smaller remaining items

- A is read three times per solve outside the ALM loop: the final
  `chi2_vector`, and two matvecs in `kkt_violation`. ~3 full passes over the
  matrix, minor now that the per-iteration pass is gone.
- `A` and `X` are both resident for the whole solve (~250 GiB at omega Cen).
  Now that chi2 comes from the residual, A is used inside `solve_adelie_alm`
  only for `A[0] @ w`. Everything else it is needed for - `chi2_vector`,
  `kkt_violation` - is expressible through `X` and `col_norm`, so A could in
  principle be released once X is built, halving the resident matrices. This
  is the largest remaining memory win and has not been attempted.

## Measured: a full adelie solve on real data (2026-08-18)

A complete `NNLS.solve()` with `nnls_solver = "adelie"` on the NGC5139
combined library (A = 322403 x 3840 = 9.2 GiB) finally ran to completion:

    orblib read   51.6s, peak RSS 15.9 GiB
    solve()     3149.9s, peak RSS 23.2 GiB
    chi2_tot=899080.802753  chi2_kin=762905.557097
    sum(w)=0.9999999998  nonzero=659 of 3840 orbits

Profile (tottime), taken BEFORE the chi2-from-residual change but WITH the
assembly and F-ordering changes:

        200  1777.7s  adelie/solver.py bvls
        200   572.1s  adelie/state.py:3272 (lambda)
          1   417.5s  weight_solvers.py solve_adelie_alm  (own time)
       1035   305.1s  numpy.ufunc.reduce
          2    67.3s  orblib.py read_orbit_base (cumulative)
          1     2.0s  weight_solvers.py construct_nnls_matrix_and_rhs

Two things worth taking from this.

**The chi2 line was ~20-25% of the whole solve.** `solve_adelie_alm`'s own
tottime of 417.5s is dominated by the `A @ w` matmul in the old chi2 line
(numpy C calls are attributed to the calling Python frame), and much of the
305s of ufunc reduces is the `np.sum` of the squared residual next to it.
The standalone microbenchmark puts the matmul alone at 2618ms x 200 = 524s.
So the residual shortcut removes roughly a fifth to a quarter of total solve
wall time here - measured, not extrapolated.

**Peak RSS is 23.2 GiB for a 9.2 GiB matrix**, which is A and X resident
together plus change. That is the direct evidence for the "release A once X is
built" item above: it would take peak from ~2.5x the matrix to ~1.5x.

Assembly is now 2.0s of a 3150s solve, i.e. no longer worth optimising.
Reproduce with `dev_tests/_real_solve_profile.py adelie`; note it takes ~55
minutes on a laptop.

## Bulk-reading the PM histograms needs no format change (2026-08-18)

Checked directly, by walking the raw record markers of a decompressed
`orblib_pm_hist.dat` (232 MiB):

    file 231.9 MiB, header record 16 bytes
      rec 0: hdr=16B ints=[8, 8, -9, -9] EMPTY (no payload)
      ...
      rec 7: hdr=16B ints=[0, -2, 0, 2] markers_match=True
             payload=5 values, expected 5, match=True

The 2D file has exactly the same shape as the 1D one that
`_walk_losvd_records` already parses, with a four-int header instead of two::

    1D:  record: [int32 ivmin][int32 ivmax]
         record: [float64 x nv]                       if ivmin <= ivmax
    2D:  record: [int32 ivmin0][int32 ivmin1][int32 ivmax0][int32 ivmax1]
         record: [float64 x nv0*nv1], Fortran order   if both ranges non-empty

So bulk-reading it is a parser change, not a format change: the same bytes,
read once into a buffer and scattered with numpy, instead of two
`scipy.io.FortranFile` record reads per (orbit, aperture) pair.

Three things make this more tractable than it looks:

- **1D and 2D histograms are in separate files** (`orblib_in` and
  `orblib_in_pm`), so the two parsers never interleave. The existing
  `all(i == 1 for i in hist_dim)` gate is stricter than it needs to be: a
  library with proper motions could already bulk-read its 1D sets today.
  That alone is small here though - the PM set has 1405 apertures against the
  LOSVD set's 163, so the PM file is where the calls are.
- `_gather_float64` and `_value_batches` are dimension-agnostic and can be
  reused unchanged. The new work is the walker stride (6 int32 slots per
  header record instead of 4) and a 2D scatter, where a payload index `t`
  maps to `(t % nv0, t // nv0)` because the payload is Fortran-ordered.
- The pairs are very sparse - 11 of the first 12 are empty, and the 1D case
  measured 77% empty - so the loop currently spends most of its time reading
  sentinels. The 24x won on the 1D path is a reasonable floor for what this
  would win, not a ceiling.

The generalisation that is actually needed: `_read_losvd_hist_vectorised`
assumes `n_pairs = norb * n_ap_total`, i.e. that every aperture is in its
file. For a mixed library each parser must be given only the apertures of its
own dimension, plus the mapping back to global aperture indices.

Correctness is cheaply checkable: run both paths on the same real library and
compare the filled histograms with `np.array_equal`.
