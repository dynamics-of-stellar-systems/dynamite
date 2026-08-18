# Weight-solve and orbit-library read performance

Work on `feature/orblib-performance`, 2026-08-18. Every number here was
measured on a real orbit library, not estimated, and every change to numerical
code was checked to produce bit-identical output.

## Setup and how to reproduce

`NGC5139_veldist_combined_bigger` (MUSE BayesLOSVD + HST proper motions, 3840
orbits, 163 LOSVD + 1405 PM apertures) in `~/research/omegacen/dynamite_dataprep`.
The assembled NNLS matrix is 322403 x 3840 float64 = 9.2 GiB. Production omega
Cen is ~45000 orbits and ~371212 constraints, ~125 GiB, i.e. about 13x this.

Run from the dataprep directory with the `main` conda python:

| script | what it does |
|---|---|
| `dev_tests/_real_orblib_check.py <out.npz>` | profiles matrix assembly, dumps A/b for revision-to-revision `np.array_equal` |
| `dev_tests/_real_orblib_read_profile.py` | profiles `read_vel_histograms()` |
| `dev_tests/_real_hist_read_check.py` | sha256 fingerprints the histograms a real library reads |
| `dev_tests/_real_solve_profile.py adelie` | profiles a full `solve()` (~55 min on a laptop) |
| `dev_tests/_real_alm_chi2_check.py <A.npz> [rows]` | runs the real ALM on a row-subsampled matrix, rechecks chi2 |
| `dev_tests/_real_chi2_kinmap_check.py` | checks the orblib comes back unmutated from assembly |

## Where a model's time goes now

A full adelie `solve()` on this library: **3149.9s, peak RSS 23.2 GiB**,
`sum(w) = 0.9999999998`, 659 of 3840 orbits nonzero. Profiled before the ALM
chi2 change but after the assembly work:

```
 200  1777.7s  adelie bvls
 200   572.1s  adelie state lambda
   1   417.5s  solve_adelie_alm (own time)   <- removed by the chi2 change
1035   305.1s  numpy.ufunc.reduce            <- mostly the same chi2 line
   2    67.3s  read_orbit_base               <- now ~13s
   1     2.0s  construct_nnls_matrix_and_rhs <- was 17.2s
```

The solver itself is now the cost, which is where it should be.

## What changed

### 1. Three redundant full-size copies of A (`a69c44e`)

- `A_normalized`/`b_normalized` were built unconditionally in `solve()` but are
  unread on the adelie path, holding a second full copy of A for the whole run.
  Now built only for the scipy/cvxopt branches, the only consumers.
- `A.astype(np.float64)` in the ALM loop is a no-op cast at the default
  `nnls_dtype`, but `.astype` copies regardless - a full copy per iteration.
  Now `copy=False`.
- `X` was `np.vstack(...)` -> `X / col_norm` -> `np.asfortranarray(...)`: three
  chained full-size allocations, peaking at ~4x the matrix with A alive. Now
  fills one F-ordered buffer in place.

### 2. Matrix assembly, 17.2s -> 3.4s, bit-identical (`860fe56`)

- **Pre-size `con`/`econ`/`orbmat`** instead of growing `orbmat` by `np.vstack`
  per kinematic set. Row counts come from
  `get_observed_values_and_uncertainties()`, which never touches the orbit
  library, so a cheap pre-pass sizes everything; its results are reused in the
  loop rather than recomputed. Peak drops from ~2.2x the matrix to ~1.4x.
- **Allocate `orbmat` Fortran-ordered.** Blocks arrive as
  `(n_orbs, n_constraints)` and are written in transposed, which is a memcpy
  into an F-ordered destination and a strided shuffle into a C-ordered one:

  | per 2.3 GiB block | C-ordered | F-ordered |
  |---|---|---|
  | block write | 1.00 s | 0.08 s |
  | divide by `econ` | 0.12 s | 0.08 s |
  | build adelie's F-contiguous `X` | 2.38 s | 0.33 s |
  | `A @ w` (chi2) | 0.04 s | 0.07 s |

- **Reshape the destination, not the source.** `orb_kins` is a `moveaxis` view,
  so `np.reshape(orb_kins, (n_orbs, -1))` cannot be a view and copied the whole
  block: 1.33s and a 2314 MiB peak, against 1.03s and **0 MiB** reshaping the
  destination instead. At omega Cen scale that temporary is ~100 GiB.

  `ndarray.reshape` silently returns a **copy** when it cannot return a view,
  which as an assignment target is a silent no-op leaving the block at zero.
  The `np.shares_memory` assert exists to make that loud; do not remove it.
  (`.shape =` raises natively but is deprecated in numpy 2.5, and the package
  supports numpy>=1.26.)

Also asserts that the orbit library and the kinematic data agree on the
constraint count per set - previously an unstated invariant.

### 3. LOSVD rebin through BLAS (`7a5e340`)

`BayesLOSVD.rebin_orblib_to_observations` contracted with
`np.einsum(..., optimize=False)`, which runs a naive loop instead of BLAS. On a
549 MiB stand-in:

```
optimize=False   2401 ms, peak alloc 220 MiB (the output alone)
optimize=True     363 ms, peak alloc 769 MiB
chunked gemm      371 ms, peak alloc 290 MiB   <- what was implemented
```

`optimize=True` is a trap: its cost model counts flops, not bytes, and it
reaches BLAS only by transposing a full copy of the largest array in the run.
Chunking over orbits gets the same speed with a bounded temporary.

`GaussHermite.get_gh_expansion_coefficients` also got `optimize=True` - only
~1.5x, but its contraction path was checked at production shapes first
(largest intermediate ~108 MB, no orbit axis).

### 4. ALM chi2 from the residual (`6fe911c`)

Every multiplier update recomputed `chi2 = sum((A @ w - b)**2)`, a full pass
over the matrix, up to `adelie_alm_iters` (default 200) times.

`X` is `A[1:]` with unit-L2 column scaling and `w = beta / col_norm`, so
`X[1:] @ beta` is exactly `A[1:] @ w`, and adelie returns
`resid = y - X @ beta`. Only row 0 of A is missing, which `X` replaces with the
ALM penalty row:

    chi2 = (A[0] @ w - b[0])**2 + resid[1:] @ resid[1:]

Measured on the real matrix: **2618 ms -> 0.04 ms** per iteration, ~520s per
model here, and 20-25% of total solve wall time in the full profile above.

`state.resid` is `y - X @ beta` unweighted (verified to 2e-15, uniform and
non-uniform weights); `state.loss` is `0.5 * sum(weights * resid**2)` and would
have silently folded in the observation weights and the penalty row.

Verified three ways: per-iteration against `A @ w` in a unit test (0 to 3e-15
relative), and end to end on real subsampled matrices -

```
rows    reported (from resid)   direct A@w - b      final gap
8000    16910.98                16910.9799656940    3.5e-09
40000   91786.1095              91786.109536        7.3e-08
```

### 5. Bulk-read the 2d proper-motion histograms (`37e4891`, `505ea6b`)

`_read_losvd_hist_vectorised`, the bulk parser behind the earlier 24x LOSVD
read speedup, was gated on the library holding nothing but 1d histograms. One
proper-motion set puts a 2 in `hist_dim` and the whole library - 1d sets
included - fell back to the per-record loop. That is exactly the omega Cen
production setup, and its PM set has 1405 apertures against the LOSVD set's
163, so nearly all of the 5.13 million record reads came from the path with no
bulk parser.

Added `_walk_pm_records` and `_read_pm_hist_vectorised`. **No file format
change**: the pm file has the same Fortran framing with a four-int header
(`[ivmin0][ivmin1][ivmax0][ivmax1]`, 6 int32 slots) and a Fortran-ordered
payload, so payload index `t` maps to `(t % nv0, t // nv0)`. `_gather_float64`
and `_value_batches` are reused unchanged. Since 1d and 2d histograms live in
separate files, both parsers now take the aperture subset for their own file
rather than assuming every aperture appears in it, and the gate becomes just
`not legacy_file and not return_intrinsic_moments`.

```
read_vel_histograms   50.9s -> 13.7s        (3.7x)
under cProfile       163.2s -> 12.5s
function calls       119,106,779 -> 198,708
```

Bit-identical: sha256 of both filled histograms matches the per-record loop.

3.7x rather than 24x because what remains is not parsing - **7.3s of the 12.5s
is `BufferedReader.read` inside bz2 decompression**. The read is now
decompression-bound, so further gains here are a storage-format question, not a
parsing one. (The 13x under cProfile overstates the change: profiling 119
million calls exaggerated the old path.)

**Trap:** the pm_hist file has NO header record, unlike losvd_hist.
`read_orbit_base` consumes a header from the losvd file only. Reusing the 1d
walker's header skip silently ate the first data record; every later record
still parsed self-consistently, one record out of phase, and it only surfaced
1.8 million pairs later as the truncation error.

## Testing

Unit tests, runnable directly, no framework:

| test | covers |
|---|---|
| `test_hist_bulk_read.py` | both bulk parsers vs a transcription of the per-record loop |
| `test_nnls_matrix_assembly.py` | pre-sized F-ordered assembly vs the old vstack version |
| `test_losvd_rebin_chunking.py` | chunked BLAS rebin vs the einsum it replaced |
| `test_alm_chi2_from_resid.py` | ALM chi2 identity, driving the real adelie bvls |

`test_hist_bulk_read.py` is the important one, because a binary parser fails
silently. 18 combinations: six layouts (1d only multi-set, 2d only multi-set,
mixed, 2d-before-1d, single aperture, wide 2d) x three conditions (half empty,
nothing empty, entirely empty), each at batch sizes 1, 7 and 1e9. The mixed
layout interleaves 1d, 2d, 1d, 2d - the omega Cen shape - so the per-file
aperture subsets are not contiguous blocks. A truncated file must raise.

`dev_tests/_mutate_hist_bulk_read.py` checks the suite can actually fail: it
breaks each load-bearing line in turn and confirms a red result. All seven
caught - walker stride, the bogus header skip, C-order payload decode, row/col
swap, dropped `idx_ap_reset`, wrong fast-axis length, and the 1d parser
ignoring its aperture subset. Two of those are mistakes actually made writing
this code.

### Residual risk

- The synthetic tests use 1-5 orbit libraries: they exercise logic, not scale.
  `_value_batches` at real sizes is covered only by the real-library check.
- Only one real library was checked, and its two halves share a layout.
- `return_intrinsic_moments=True` and the legacy format still use the old loop.
  Unaffected, but also unverified here.
- Nothing has been run at full omega Cen scale. The highest-value next check is
  `_real_hist_read_check.py` on a production library, diffed against a
  pre-change revision.

## Considered and deliberately not done

- **`Histogram2D.rebin_orblib_to_observations`** (the PM rebin) looks terrible -
  nested Python loops recomputing `overlapping_bins` inside the outer loop -
  and the overlap weights are exactly separable, so it collapses to two
  tensordots. Written and verified to 3e-16, but only 1.0-1.7x faster at 20x20
  to 60x60 velocity bins, because the overlaps stay local so the loop is
  effectively linear in bin count. In the test config it does not even run
  (`No rebinning necessary`). Revisit only for much finer PM grids.
- **`chi2_kinmap` reusing the orblib.** `chi2_kinmap` re-reads the whole library
  because assembly zeroes the edge velocity bins in place; restoring those two
  slices avoids the re-read (106.5s on this library). Implemented and verified
  (the orblib comes back `np.array_equal`), then dropped from the branch: it
  only helps all-GaussHermite configurations, since `chi2_kinmap` returns nan
  before reading anything if any set is not GaussHermite, so it does nothing
  for omega Cen while still touching the hot assembly path. Preserved on branch
  `chi2-kinmap-orblib-reuse` for upstream.

## Open

- **Release A once X is built.** Peak RSS is 23.2 GiB for a 9.2 GiB matrix -
  A and X resident together. After the chi2 change, A is used inside
  `solve_adelie_alm` only for `A[0] @ w`. Everything else it serves is
  expressible through `X` and `col_norm`: `A[1:] @ w = X[1:] @ beta`,
  `A^T r = a0*r0 + col_norm * (X[1:]^T r_rest)`, and
  `||A_col||^2 = col_norm^2 - mu + a0^2`. Two obstacles: A is a local in
  `solve()` so the callee cannot release it, and `kkt_violation` is the only
  optimality certificate in the code, so rewriting it against a scaled
  surrogate deserves care. Would take peak from ~2.5x the matrix to ~1.5x.
- **A is still read three times per solve** outside the ALM loop: the final
  `chi2_vector` and two matvecs in `kkt_violation`. Minor now.
- **`nnls_solver: "cvxopt"`** builds `P = A_normalized.T @ A_normalized`, an
  (n_orbs, n_orbs) Gram matrix. At 45000 orbits that is 16 GB and ~7.5e14
  flops. Unusable at production scale; fine for small problems.
- The configs here set `VECLIB_MAXIMUM_THREADS=1` to dodge an Accelerate
  SIGSEGV in `scipy.optimize.nnls`, which also serialises every BLAS call in
  the solve. Worth remembering when timing the scipy path.

## Cherry-picking

Checked by applying each commit's patch to `master` in a scratch worktree:

```
CLEAN     perf(weight-solvers): stop making three redundant full-size copies of A
CLEAN     perf(weight-solvers): pre-size and F-order the NNLS matrix
CLEAN     perf(kinematics): route the LOSVD rebin through BLAS
CONFLICT  perf(weight-solvers): take ALM chi2 from adelie's residual
```

The conflict is an ordering artifact, not a logical dependency: the ALM chi2
change replaces the line the first commit last touched. The identity holds
regardless; applying it alone needs a trivial textual fixup.

## Corrections made along the way

Recorded because the wrong versions were stated confidently before being
measured:

- The hour-long assembly phase was blamed first on `np.vstack`, then on
  `transform_orblib_to_observables`. Neither. The transforms are noise (LOSVD
  rebin 0.13s, PM transform 0.002s) and the cost was in plain assignment
  statements, i.e. memory layout. Scaling 17s by 13x puts omega Cen assembly at
  3-4 minutes even before these changes, so whatever took an hour was not this
  function.
- Local solve runs were reported as having crashed. They had not - they were
  killed by how they were launched, and they complete.
