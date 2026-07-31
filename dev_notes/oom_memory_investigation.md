# NGC5139 OOM investigation: findings, what shipped, what's planned

Investigation triggered by a production OOM (peak RSS 964.63 GB at
12:53:34, 18 concurrent weight-solving workers at 44-68 GB each, against
1416.82 GB total system RAM). Source data: `~/Downloads/dynamite.log` and
`dynamite_monitor.jsonl` from that run (not checked in).

## Root cause

Not a duplication bug: log-verified exactly one `NNLS/adelie` solve per
model (27 solves = 9 orbit libraries x 3 `ml` values, one orbit library
reused across `ml` scalings without re-integration).

The per-worker memory driver is `weight_solvers.py`'s `NNLS.solve()`:
`orblib.read_vel_histograms()` sets `orblib.vel_histograms`/
`intrinsic_masses`/`projected_masses`, and `orblib` stays referenced (kept
alive by the caller) through the entire ALM solve, even though
`construct_nnls_matrix_and_rhs` fully copies their contents into `A`/`b`
early on and nothing touches them again afterward.

## Tried and ruled out: `del`-ing the retained orblib data

Measured directly in the real pipeline (NGC6278 test orbit libraries built
at 0.36 GB, 1.2 GB, 3.5 GB, and 7.7 GB): explicitly `del`-ing
`orblib.vel_histograms`/`intrinsic_masses`/`projected_masses` right after
`construct_nnls_matrix_and_rhs` reclaimed **~0% of RSS at every size**,
confirmed not to be a Python reference leak (`gc.get_referrers` showed only
the expected `Histogram` wrapper). Root cause: heap fragmentation from
*allocation order*, not size - other same-size temporaries in the very same
function call **did** get reclaimed cleanly, ruling out a pure size
threshold. Two independent research passes confirmed this is a known NumPy
anti-pattern (building large arrays through generations of `vstack`/
`concatenate`/mirror-copy rather than pre-allocating and writing in place),
and that it is not macOS-specific - glibc on Linux has the same
"hole in the middle of the heap" failure mode for freed memory that isn't
at the top of the break or in its own mmap'd region.

Also ruled out: sparse matrix representation for `orbmat` - the matrix
itself is only ~2.76 GB against a 44-68 GB per-worker budget; not worth a
much larger, riskier rewrite (dense allocation happens too early in
`orblib.py`/`kinematics.py` to exploit cheaply).

## Shipped

1. **`ncpus_weights` resized** to a value sized to worst-case per-worker RSS
   against a 50%-of-total-RAM budget (already done before this note was
   written).
2. **`nnls_dtype: 'float32'` option** (`weight_solver_settings`, optional,
   default `'float64'`) - `dynamite/weight_solvers.py`. Downcasts
   `orblib.vel_histograms`/`intrinsic_masses`/`projected_masses` right after
   `read_vel_histograms()`, and builds `con`/`econ`/`orbmat`/the ALM solve
   arrays in the configured dtype. Validated three ways on real NGC6278
   orbit libraries (1440-orbit original and a 4800-orbit scratch build):
   chi2 agrees with float64 to <0.001% in every test, weight-vector
   correlation >=0.99997, KKT violation stays within the range of
   previously-accepted float64 solutions (never close to the 0.1 warning
   threshold). Measured RSS reduction: matrix-only ~24%, full end-to-end
   (orblib downcast + matrix + solve, two fresh processes compared directly)
   ~24-28% on the tested sizes. Not yet validated on a dataset whose
   constraint-row scaling differs from NGC6278/omega Cen - `adelie_mu` is
   documented to have the same caveat.

   Required an upstream fix in the process: adelie's `bvls()` crashes on
   float32 input via its default `weights=None` path (`weights` defaults to
   `np.full(n, 1/n)`, always float64, and the internal state constructor's
   strict `copy=False` downcast then raises under numpy>=2.0). Worked around
   by passing an explicit float32 `weights` array. Separately hit a NEP 50
   gotcha: `np.sqrt(mu)` returns a numpy float64 *scalar* (not a Python
   float), and multiplying a numpy scalar against a float32 array forces
   upcast under NumPy 2.0's promotion rules (a plain Python float literal
   would not) - fixed via `np.full(..., sqrt_mu, dtype=dtype)` instead of
   `sqrt_mu * np.ones(...)`.

3. **`ncpus_weights_maxtasksperchild` option** (`multiprocessing_settings`,
   optional, default `None` = previous behaviour) -
   `dynamite/model_iterator.py`'s `SplitModelIterator`/`ModelInnerIterator`.
   Recycles each weight-solving `Pool` worker after N models instead of
   reusing it for the whole run. This guards against RSS creep *across*
   models in a long-lived worker (since process exit unconditionally
   returns everything to the OS regardless of fragmentation state) - it does
   **not** reduce memory used *during* one model's solve, which was the
   original ask; kept as a complementary hygiene fix, not the primary lever.

## Also shipped: orbit-library copy-chain preallocation

See `docs/superpowers/plans/2026-07-31-orblib-copy-chain-preallocation.md`
for the design. `duplicate_flip_and_interlace_orblib` + `combine_orblibs`
(which built the final per-kinematic-dataset array through 3 full-copy
generations, ~3-4x peak) were replaced with `combine_and_mirror_orblibs`, a
single function that allocates the final array once and writes tube/box
data directly into its slices. Prototyped and validated first: bit-exact
equivalence, full weight-solve output bit-identical (A/b/weights/chi2),
dtype-preserving, ~4% faster, peak RSS reduced 4.7-12.1% (fresh-process
comparisons, small and 1.2GB test orbit libraries). Then landed in
`orblib.py` and re-verified against the real patched file (three
independent scripts, exact pre-change chi2 reproduced). Composes with
`nnls_dtype` (different axis: allocation count vs. bytes per element) and
sidesteps the reclaim question entirely since it reduces what's allocated,
not what's freed.

Not yet tested against a real orbit library: `mirror=False` (bar/disk
systems) and 2D histograms (proper motions) - both branches are
implemented but no locally-built orbit library exercises either case.

## Found (and fixed, carefully) a real double-read in chi2_kinmap

Running the "run all tests" pass at full scale surfaced something new:
`NNLS.solve()`'s `chi2_kinmap()` step (run whenever all kinematics are
GaussHermite) calls `Analysis.get_gh_model_kinematic_maps()` once per
kinematic set, and that method does `orblib = model.get_orblib()` - and
`Model.get_orblib()` (`model.py:1184`) always constructs a brand-new
`LegacyOrbitLibrary` with no caching. So a full, independent
`read_vel_histograms()` happens once per kinematic set, on top of the read
already done for the main NNLS solve - N+1 total reads for N kinematic sets
per model, not 1.

The obvious fix - pass the caller's already-read `orblib` into
`chi2_kinmap`/`get_gh_model_kinematic_maps` instead of fetching fresh ones -
is **unsafe as a naive reuse** and was caught by testing before it shipped:
`construct_nnls_matrix_and_rhs` (`weight_solvers.py:790-791`) zeroes the
first and last velocity bin of each 1D histogram *in place* on
`orblib.vel_histograms`, unconditionally (not gated by the `CRcut` setting),
as a side effect of building the NNLS matrix. Before any fix, `chi2_kinmap`
happened to get an unmutated orblib (a fresh one, never passed through
`construct_nnls_matrix_and_rhs`). Reusing the *outer*, already-mutated
orblib changed `chi2_kinmap`'s value silently - confirmed by A/B testing
with `git stash`: chi2_kinmap went from 695753.906576 (correct, matching
the pre-change baseline) to 650393.493585 (wrong) on the same weights/
chi2_tot/chi2_kin, on the small NGC6278 test model.

**The shipped fix**: `chi2_kinmap` fetches its own fresh, unmutated orblib
**once** (not per kinematic set) and reuses that single fresh read across
the loop over kinematic sets - deduplicating N reads down to 1, without
ever touching the outer (mutated) orblib. Verified this reproduces the
exact pre-change chi2_kinmap value.

**Validated on a real N>1 case**: the NGC5139 (omega Cen) "mini" orbit
library at `~/research/dynamite_analysis/test_head` (already integrated,
3 GaussHermite kinematic datasets, no `all_models.ecsv` existed yet so one
was built by registering the fixed parset directly) gives a genuine N=3
test. `git stash`-based A/B comparison:

| | pre-fix | with fix |
|---|---|---|
| chi2_tot | 1739592.046419 | 1739592.046419 |
| chi2_kin | 118613.274617 | 118613.274617 |
| chi2_kinmap | 133142.838583 | 133142.838583 (exact match) |
| solve() wall time | 7.04 s | 5.81-5.90 s (~17% less) |

Confirms both correctness (bit-identical chi2_kinmap with N=3, not just
N=1) and a real, measurable time saving from eliminating 2 of the 3
redundant per-kinematic-set reads. This "mini" library is still small; on
the full production omega Cen orbit library (~34.7 GB combined histograms
per the production memory note) each eliminated read would cost far more
in both wall-time and transient peak memory, so the production benefit
should be substantially larger than this 17%. `float32` (`nnls_dtype`) also
verified working correctly on this same N=3 model: chi2_tot=1739592.021756
vs the float64 value 1739592.046419, agreement to <0.00001%.

## Caveats

All memory measurements were on macOS. The pre-existing production memory
note (see `dynamite_production_memory` in agent memory) already flagged
that a prior macOS extrapolation overestimated a comparable number by 2x on
the real Linux cluster - re-measure there before trusting these exact
percentages, though the *mechanism* (fragmentation from allocation order)
is expected to generalize, not to be macOS-specific.
