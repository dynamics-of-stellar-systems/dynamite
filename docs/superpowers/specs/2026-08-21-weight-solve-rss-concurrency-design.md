# Weight-solve resident memory: fused-X assembly, streamed reads, concurrency — Design

Date: 2026-08-21
Branch: `feature/orblib-performance` (continues the 2026-08-18 sprint)
Target machine: 192-core / 1416 GB KVM node (`PM_grid` cluster), no local scratch, no GPUs.

## Motivation

The Aug-18 sprint made every stage of an omega-Cen weight solve faster
(assembly 17.2s -> 3.4s, ALM chi2 2618ms -> 0.04ms/iter, histogram reads
50.9s -> 13.7s) and removed three redundant copies of A, but the *resting*
resident set of a production solve is still dominated by two irreducible-looking
objects held simultaneously:

- `A`, the error-normalized NNLS matrix (~371212 x ~45000 float64, ~120 GiB),
- `X`, adelie's augmented column-scaled copy of A[1:] plus the penalty row
  (another ~120 GiB),

plus all kinematic sets' histograms (~35 GB class, extrapolated, never measured
here) held from `read_vel_histograms()` through assembly. Peak is therefore of
order 2.5x the matrix; at production scale that caps the node at roughly 3-4
concurrent weight solves before OOM risk, and `ncpus_weights: 1` is currently
the safe setting. The production grid wants as many concurrent solves as fit.

The macOS OOM investigation found `del` reclaimed ~0% after assembly, but that
was never re-tested on Linux glibc, where huge numpy allocations are mmap-backed
and usually return to the OS on free.

## Goal

Minimize per-worker peak RSS during a production-shape adelie weight solve so
that as many concurrent solves as possible fit on the 1.4 TB node, with these
constraints agreed up front:

- Algorithm changes are in scope; correctness gates are exact reproduction at
  float64 (bit-identical weights/chi2 against the recorded baseline) and,
  for float32, the already-demonstrated <0.001% chi2 tolerance.
- scipy/cvxopt solver paths keep their existing code path untouched.
- No out-of-core tiering: the node has no usable scratch (25 GB root disk),
  so everything must fit in RAM.

## Approach (approved)

Profile-led structural sprint in three phases, each gated by measurement:

1. Baseline RSS profile of one real production-shape solve on this node.
2. Two structural changes: fused-X assembly (never materialize A) and
   streamed per-set histogram reads.
3. Stack the validated float32 mode; re-profile; publish the concurrency table.

An out-of-core / reformulated-ALM approach was considered and parked: it has
the highest ceiling but abandons the validated adelie codepath for random
column access over network FS, and X is already near-irreducible without it.
Config-only tuning (float32 + recycling alone) was rejected as insufficient
(~25-30% cut, ~2 concurrent solves).

## Phase 1 — Baseline profiling campaign

New script `dev_tests/_real_solve_rss_profile.py` (modeled on the existing
`_real_solve_profile.py`):

- In-process phase markers around: orblib read (per kinematic set), matrix
  assembly (per block), X build, ALM iterations (logged every 10th),
  kkt_violation, final chi2_vector.
- A 1 Hz sampler thread reading `/proc/self/status` (VmRSS, VmHWM) and
  `/proc/self/smaps_rollup` (Anonymous, Pss), writing CSV compatible with the
  PM_grid monitor CSVs.
- A standalone reclaim probe: allocate histogram-sized and matrix-sized
  arrays, free them, sample RSS after `del` + `gc.collect()`, optionally
  `malloc_trim(0)` via ctypes; report reclaimed percentage on Linux.

Run configuration: xeast input pack, orblib_000_000 hardlinked into a scratch
output directory, fixed parset ml02.60, adelie solver, float64 - identical to
the completed baseline run whose chi2_tot=2770835.03357815 /
kinchi2=335126.55535470234 is recorded in
`PM_grid/NGC5139_adelie_xeast_output/all_models.ecsv`. That number doubles as
the end-to-end ground truth for Phase 2 validation.

Two modes: `--alm-iters N` capped run (~30 iterations, about 1 hour; per-
iteration memory is flat so the plateau is visible early) for iteration speed,
and a full uncapped run for the final record.

Questions the profile must answer before Phase 2 proceeds:

1. Measured bytes of each kinematic set's histograms and the masses arrays at
   this scale (the ~34.7 GB figure was extrapolated).
2. Does `del` reclaim on Linux glibc here, post single-allocation combine?
   If yes, streamed frees are reliable; if no, Phase 2's streaming still helps
   (it avoids co-residency rather than relying on reclaim) but float32 becomes
   more important.
3. Where peak sits now: read/assembly overlap vs X build vs ALM plateau.
4. Transient during `combine_and_mirror_orblibs` at production scale.

Deliverable: `dev_notes/weight_solve_rss_profile.md` with the phase-by-phase
GiB table, timeline CSV(s), reclaim results, and the measured per-worker budget.

## Phase 2a — Fused-X assembly (never materialize A)

Today `construct_nnls_matrix_and_rhs` builds A, then `solve_adelie_alm`
re-copies A[1:] into an F-ordered augmented X and column-scales it in place -
two full matrices resident together for hours.

Change: assemble the augmented matrix directly.

- Allocate `X = np.empty((n_rows, n_orbs), dtype, order='F')` once - exactly
  `n_rows` rows, because the `sqrt_mu` penalty row *replaces* the total-mass
  row rather than being added to it (`X[0, :] = sqrt_mu`).
- Stream each block (total-mass row, intrinsic masses, projected masses, one
  block per kinematic set) straight into `X[1:, :]`, applying its `econ`
  division as written - same per-element operation order as today's post-hoc
  divide.
- Compute `col_norm = np.linalg.norm(X, axis=0)` afterwards exactly as today,
  then divide X in place. Deliberately NOT an incremental per-block norm
  accumulation: that would change floating-point rounding. This way X,
  col_norm and y are bit-identical to the current two-step construction by
  construction.

`solve_adelie_alm` is refactored to accept `(X, col_norm, b_rest, row0)`
instead of `(A, b)`; A ceases to exist on the adelie path. Consumers of A are
rewritten against identities already derived in
`dev_notes/weight_solve_performance.md` ("Open"):

| consumer | today | after |
|---|---|---|
| ALM chi2 | already from state.resid + row-0 dot | unchanged |
| chi2_tot / chi2_kin | `(A@w - b)**2` | row-0 term + `resid[1:]**2`; identical vector because `resid[1:] == b_rest - A_rest@w` in the error-normalized space chi2 uses |
| kkt_violation(A, b, w) | two matvecs over A | surrogate via `A^T r = a0*r0 + col_norm*(X[1:]^T r_rest)` etc.; unit-tested against direct computation |

Peak during solve drops from A+X (~2.4x matrix) to X alone (~1.2x). The
scipy/cvxopt paths keep calling the existing A-building function unchanged.

Guards: keep the `np.shares_memory` assert pattern on every block write
(ndarray.reshape silently returning a copy would leave blocks at zero); keep
the zero-column col_norm guard.

## Phase 2b — Streamed per-set histogram reads

Today `read_vel_histograms()` materializes all kinematic sets' histograms up
front and holds them through assembly. Change: assembly drives reads set by
set - read set i's histograms, transform to observables, write the block into
X, free, then read set i+1. Only one set's histograms plus accumulated-X are
ever co-resident.

- Requires a per-set read API on the orbit library. The bulk parsers already
  accept per-file aperture subsets, so this is a thin refactor, not new
  parsing.
- The in-place edge-velocity-bin zeroing done during assembly becomes
  irrelevant on this path (histograms die inside the loop); `chi2_kinmap`'s
  fresh-read contract is untouched and still returns nan before reading
  anything when any kinematic set is not GaussHermite (the omega Cen case).
- Intrinsic/projected masses stay resident (small).

If Phase 1 shows Linux reclaims cleanly, streaming mainly shortens the
co-residency window; if not, streaming is what makes the win real.

## Phase 3 — float32 stack, validation, deliverables

`nnls_dtype=float32` composes with both changes (dtype plumbing, the adelie
`weights=` kwarg workaround and the NEP 50 scalar gotcha are already shipped
and validated). Expected landing zone pending Phase 1 numbers: ~60 GiB
per worker at float32 fused-X, i.e. roughly 8-12 concurrent solves on this
node versus 3-4 today.

Validation gates, each change independently:

1. Unit: fused assembly vs the existing two-step on synthetic shapes -
   bit-identical X, col_norm, y.
2. Unit: surrogate chi2/kkt identities vs direct computation on row-subsampled
   real matrices (extending `_real_alm_chi2_check.py`'s method).
3. Real library: extend `_real_orblib_check.py` to dump X/col_norm/y for
   revision-to-revision `np.array_equal` comparison.
4. End to end at float64: rerun the xeast ml02.60 solve; must reproduce
   chi2_tot 2770835.03357815, kinchi2 335126.55535470234, and weights exactly.
5. Re-profile after each change with the Phase 1 tooling; final deliverable is
   a concurrency table (measured per-worker peak x {float64, float32} x
   {streamed}) with the recommended `ncpus_weights`.

## Testing

New dev_tests following existing conventions (direct-run scripts plus pytest-
collectable modules):

- `test_fused_assembly.py` - gate 1 across shapes incl. empty sentinels,
  mixed 1d/2d layouts, batch extremes.
- `test_surrogate_chi2_kkt.py` - gate 2.
- Extended `_real_orblib_check.py` / new fingerprint script - gate 3.
- Reclaim probe results recorded in the notes regardless of outcome.

Notes consolidated into `dev_notes/weight_solve_performance.md`; the RSS
profile lives in `dev_notes/weight_solve_rss_profile.md`.

## Risks and open questions

- kkt_violation is the only optimality certificate in the codebase; the
  surrogate rewrite is unit-tested against direct computation and its scaled
  value must remain comparable (same [0,1] definition, same reference table).
- adelie's bvls receives X directly today; the refactor keeps the same call
  signature semantics (F-contiguous, column-scaled, warm-startable states), so
  solver behaviour should be bit-identical given identical inputs - asserted in
  gate 1/4 rather than assumed.
- Fragmentation behaviour differs across libc versions; the Phase 1 probe is
  re-run after each structural change so the concurrency table reflects
  reality, not extrapolation.
- The capped-iters profiling mode changes the chi2 trajectory (fewer
  iterations); it is used only for memory observation, never for numerical
  comparison.
