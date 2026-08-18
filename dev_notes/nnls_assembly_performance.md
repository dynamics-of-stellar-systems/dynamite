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
