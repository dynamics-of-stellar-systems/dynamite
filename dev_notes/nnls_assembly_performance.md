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

NOT verified at production scale. `dev_tests/_real_alm_chi2_check.py` runs
`solve_adelie_alm` on a row-subsampled real matrix and recomputes chi2 from
the returned weights, but every attempt to run it on this laptop died without
a traceback shortly after the matrix was built - consistent with the
Accelerate/BLAS crash the configs already work around. Run it on the cluster.

## Possible follow-up: dropping A during the solve

After this change A is used inside `solve_adelie_alm` only for `A[0] @ w`
(a single row). If `chi2_vector` and `kkt_violation` were also expressed
through `X` and `col_norm`, A could be released once X is built, halving the
resident matrix count from two to one (~125 GiB). Not attempted: it changes
what `kkt_violation` is handed, and that function is the only optimality
certificate in the code.
