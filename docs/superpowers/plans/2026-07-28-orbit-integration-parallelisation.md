# Parallelising orbit integration by splitting the orbit range

Status: **feasibility tests executed, all pass**. Created 2026-07-28.
Results in section 8; the plan below is kept as written for context, and where
a test refuted an expectation that is called out rather than edited away.

Goal: cut the wall-clock cost of building an orbit library by running N
processes over disjoint orbit ranges, without touching the integrator physics.

---

## 1. Why this, and why now

Profiling a real 9-day, 130-model omega Cen grid (`dynamite_monitor.jsonl`)
gives the following split of the work actually done:

| Stage | CPU-hours | Share |
|---|---|---|
| Python (driver + NNLS) | 2039 | 57% |
| Fortran orbit integration | 1519 | 43% |
| `bunzip2` | 11.0 | 0.3% |
| `bzip2` | 3.3 | 0.09% |

The Python bucket was dominated by the scipy NNLS solve. The adelie/ALM
solver (`feature/adelie-alm-solver`) removes roughly a factor of 10 from it,
and the vectorised LOSVD read (first commit on this branch) removes a
further ~230 s per model-ml. Both of those land on the Python side, so the
pipeline is now **Fortran-bound**: orbit integration is the largest remaining
line item by a wide margin.

Compression is *not* worth attacking. An early hypothesis that `bzip2` cost
~36 min per library came from reading file mtimes; measured throughput is
28.8 MB/s, so 3.84 GB takes about 2 minutes, and the monitor puts all
compression at 0.09% of the run. `zstd -3 -T0` does beat `bzip2 -9` outright
on this data (0.18 s vs 20.84 s, and *smaller*: 183 MB vs 187 MB), so it is a
free change if that code is touched anyway, but it buys no wall clock and
would invalidate every existing `.bz2` library.

## 2. Baseline

From the Fortran's own `cpu_time` instrumentation, omega Cen:

| | CPU time | orbits | mean/orbit |
|---|---|---|---|
| tube | 12,284 s (3.41 h) | 12,000 | 1.02 s |
| box | 21,122 s (5.87 h) | 12,000 | 1.76 s |

Total **9.28 CPU-h**, wall **5.87 h**. The only parallelism today is tube ‖ box
(`orblib.py:write_executable_for_all_orbits` backgrounds two shell jobs), and
they are unbalanced, so one core idles for 2.46 h.

Measured properties that make orbit-range splitting attractive:

- **Load balance is excellent.** Slowest single orbit is 0.04% of the total;
  the top 1% of orbits account for 3.1%. There are no stragglers, so scaling
  should be close to linear.
- **Per-process overhead is negligible.** A 1-orbit run costs 0.97 s wall and
  118 MB RSS (omega Cen; 0.23 s / 83 MB for NGC6278). At N=64 that is ~7.5 GB
  and 0.06% overhead against a 522 s chunk.

Projected wall clock at 9.28 CPU-h: **11 cores -> 0.84 h (7x)**, 32 -> 0.29 h
(20x), 64 -> 0.145 h (40x).

## 3. What already exists

The Fortran **already accepts an orbit range on stdin**. In `orblib.in`:

```
line  6:  1     [starting orbit]
line  7: -1     [orbits to integrate; -1 --> all orbits]
```

read by `integrator_setup` into `integrator_start` / `integrator_number`. The
output filenames (qgrid, losvd hist, orbit classification) are also per-run
input lines, so each chunk can write its own files with no collisions. There is
even a `.tmp`-based resume mechanism, so restartability was designed in.

## 4. Known blockers

### 4.1 The range feature is broken (2-line fix, verified)

`integrator_integrate` compares an **absolute orbit index** against a **count**:

```fortran
if (dith <= integrator_dithering**3 .and. &
    integrator_current <= integrator_number) then
...
if (integrator_current > integrator_number) alldone = .true.
```

With `start=241, number=240`, `integrator_current` begins at 241, immediately
exceeds 240, and zero orbits are integrated (followed by a crash deallocating
an unallocated `vel_old`). It only works for `start=1`, where "index <= count"
coincidentally holds.

Fix, verified by rebuilding and re-running: compare against
`integrator_start + integrator_number - 1` in both places. Patch is small enough
to reproduce inline; it was reverted from the tree pending this plan.

### 4.2 The RNG makes chunks non-reproducible (the real question)

Each orbit draws `rnd = ran1(1)` to randomise its sampling phase ("start
storing the orbit after 1+? steps to avoid aliasing"), and `ran1` is a **single
global sequential stream**. An orbit's offset therefore depends on how many
orbits were integrated before it.

Measured on NGC6278 (480 orbits, split 240/240):

| Comparison | Result |
|---|---|
| orbits 1-240, monolithic vs chunk starting at 1 | **bit-identical** |
| orbits 241-480, monolithic vs chunk starting at 241 | **differ**: rel L1 14.6%, per-orbit mass 0.31% |

This is sampling noise rather than physics -- 50,000 points across
152 apertures x 203 bins is ~1.6 points per cell, so large per-bin scatter is
expected, and the integrated per-orbit mass holds to 0.31%. But it is big
enough to move chi2, and it would break reproducibility against every existing
library.

The intended fix is to seed the phase draw deterministically from the orbit
index, which makes chunking exact *and* makes results independent of how the
work was divided -- a stronger reproducibility guarantee than today. Designing
that is Test 1 below and gates everything else.

## 5. Test plan

Ordered by how fast each could kill the approach. Tests 2-4 are verification;
Test 1 is the only one with genuine uncertainty.

### Test 1 -- per-orbit RNG seeding (gating)

Change the phase draw to seed from the orbit index, then:

- **1a** chunked vs monolithic, *for every chunk*, must be bit-identical.
  Pass/fail; this is the "is chunking exact" test.
- **1b** new seeding vs current seeding, monolithic, to quantify how far the
  library moves from today's references.
- **1c** chi2 from 1b against the stored NGC6278 values. Does the science
  answer move beyond sampling noise?

If 1a passes and 1c is small, the rest is engineering. If 1c is large, every
existing reference needs re-validation and the cost calculus changes.

### Test 2 -- merge correctness

`merge(chunks)` must be bit-identical to monolithic across **all three** output
streams: `_qgrid.dat`, `_losvd_hist.dat`, `.dat_orbclass.out`. Includes
rewriting the header orbit count, and asserting that concurrent chunks do not
collide over the `.tmp` resume files (they are keyed on output filename, so
distinct names should be safe -- assert it rather than assume). The record
parser added by the vectorised-read commit handles verification.

### Test 3 -- end-to-end chi2

Build an NGC6278 library at N=4, run the weight solver, compare against the
stored references (111333.37 / 13301.23 / 99728.43). This is the acceptance
test; it is the same one that validated the vectorised read.

### Test 4 -- scaling curve

Wall time vs N on NGC6278 (480 orbits) for N = 1, 2, 4, 8, 11. The load profile
predicts near-linear; measure where I/O contention bites, since that sets the
recommended N per cluster node.

### Test 5 -- compiler portability (independent of the above)

`legacy_fortran/Makefile:59` builds with `-ffast-math -O3 -march=native
-fomit-frame-pointer -funroll-loops` (`FAST=TRUE` is set at line 9, so these
are live).

- `-march=native` risks SIGILL on heterogeneous compute nodes. Rebuild with a
  portable target, confirm chi2 unchanged, measure the speed cost.
- `-ffast-math` permits FP reassociation that can shift DOP853 trajectories.
  Measure the chi2 delta rather than assuming it is zero.

Needed before any large grid regardless of whether chunking proceeds.

### Test 6 -- script-level interactions

The combined tube script also invokes `triaxmass` / `triaxmassbin`. Chunking
must run those once, not N times. Cheap to check, annoying to find late.

### Test 7 -- determinism under concurrency

Run the same chunked config twice; output must be bit-identical. Guards against
shared-file races that a sequential test would not expose.

### Final gate

Full omega Cen chunked end-to-end, compared against the existing library. Costs
9.28 CPU-h, so run once, after NGC6278 passes 1-4.

## 6. Deliberately out of scope

- **Rewriting the integrator** (JAX or otherwise). Tracked separately in
  `dev_notes/jax_orbit_integration.md`; deferred by choice.
- **OpenMP inside the Fortran.** The integrator carries `save`d state across
  calls (`dith`, and `XOUT`/`count` inside the DOP853 output callback) on top of
  module-level histogram and qgrid arrays. Threading that means auditing every
  saved and module variable for `threadprivate`, under `-ffast-math`. Separate
  processes get the same scaling with far less risk.
- **Switching compression to zstd.** Measured; not worth the library
  invalidation on its own (section 1).

## 7. Loose ends found along the way

- **Fixture staleness.** `dev_tests/NGC6278_output` was built at
  `nE:10, nI2:8, nI3:6`, but commit `f0c20b6` restored the stock config
  (`6/5/4`), so `test_nnls_adelie_vs_scipy.py` now fails at the qgrid header
  check until the two agree. The omega Cen *mini* tree
  (`NGC5139_mini_output`) is staler still: its `infil/parameters_pot.in` is
  missing the trailing `H` line the current binary expects, so the Fortran
  cannot run there at all. The full `NGC5139_output` tree is current.
- **Memory.** LOSVD histograms are dense, not sparse: ~35 GB for a combined
  omega Cen library. That bounds how many models can be read concurrently on a
  node, independently of orbit-integration parallelism.
- **Per-ml re-reads.** `ml` only rescales the velocity axis
  (`orblib.py:1241`), so the LOSVD parse is `ml`-independent yet repeated for
  every `ml`. Caching it would nearly eliminate the read, at the cost of holding
  those tens of GB across models.

---

# 8. Results

All seven tests executed on NGC6278 (480 orbits, tube + box). Harness in
`~/research/dynamite_analysis/` scratch; the merge utility is
`merge_chunks.py`.

## 8.1 Test 1 -- chunk equivalence: PASS, after a second blocker

Chunked output is bit-identical to a monolithic run for every chunk. Reaching
that needed **two** order dependencies removed, not the one this plan
predicted:

1. **The RNG** (section 4.2), as expected. Re-seeded per orbit index. This
   needed more than a plain re-seed: seeds `base + 7919k` leave the
   generator's *first* draw U-shaped over [0,1) instead of uniform (measured
   KS p=0.02 over 2000 orbits), which would have biased every orbit's sampling
   phase the same way. Sixteen warm-up draws restore uniformity (p=0.40).
2. **`real(kind=dp), save :: stepsize` in `real_integrator`** -- not
   anticipated here. DOP853's converged step was carried into the *next*
   orbit as its initial guess, so orbit k's trajectory depended on orbit k-1.
   This is exactly the hidden-state class that section 6 cites as the reason
   to avoid OpenMP; it defeats process chunking too. Lifted to module scope
   and reset per orbit.

## 8.2 The initial step matters, and the integration is not converged

Resetting the step to zero (letting DOP853 estimate it) shifted chi2 by
**+0.45% / +2.26% / +0.27%**, i.e. 5.4 and 3.3 sigma against seed scatter.
An ablation isolated the cause to the stepsize, not the RNG. Giving DOP853 a
deterministic guess of the right magnitude (the sampling interval, `RPAR(2)`)
moved chi2 the *other* way, -2.12% on the same model.

Tightening the tolerance does **not** resolve this. chi2 wanders
non-monotonically:

| tol | old (inherited step) | new (deterministic) |
|---|---|---|
| 1e-5 | 13308 | 13049 |
| 1e-6 | 14106 | 12771 |
| 1e-7 | 12848 | 13015 |
| 1e-8 | 13331 | 13462 |

Spread across tolerance and scheme is ~1335 in absolute chi2 for **all three**
models (1.2%, 10.1%, 0.8% relative), against a seed scatter of only 30-110
(0.1-0.2%). No integrator retries fire at any of these tolerances, so the
retry logic is not the mechanism.

The sensitivity is spread over essentially every orbit, not a chaotic
minority: between 1e-7 and 1e-8, only 1 of 480 tube orbits and **0 of 480**
box orbits change by less than 1%, and the ten worst orbits carry just 8-10%
of the total difference. Box orbits are systematically more sensitive than
tube (median per-orbit relative L1 0.125 vs 0.084), as expected for the more
chaotic family in a triaxial potential.

It only partly cancels in differences: across all settings the spread in
chi2 differences between models is 630 (model0, 53% cancellation) and 862
(model2, none).

**Reading.** The orbit library is a Monte Carlo object whose per-bin LOSVDs
are reshuffled by any perturbation -- seed, tolerance, initial step, and
formerly the position of an orbit in the run. Changing the initial-step scheme
perturbs it by about as much as changing the tolerance by one decade. So
removing the order dependence is not a degradation, but it is a **one-time
re-baselining**: libraries built after this change are not comparable at the
sub-percent level with libraries built before it. Model ranking was preserved
in every configuration tested.

This is a pre-existing systematic in the code, surfaced rather than introduced
by this work, and it is worth understanding on its own terms before it is used
to set confidence intervals.

## 8.3 Tests 2, 3, 7 -- merge, end-to-end chi2, determinism: PASS

- **Merge** is byte-identical to monolithic on all four streams. One trap: each
  chunk closes its file with a 1-byte Fortran record, so a merged file must
  carry exactly one of them, not one per chunk.
- **chi2** from a merged 4-chunk and 8-chunk library is identical to the
  monolithic value to the cent (111288.54 / 13076.08 / 99634.66).
- **Determinism**: two identical chunked runs, executed concurrently with other
  jobs, produced identical output.

## 8.4 Test 4 -- scaling

| chunks | processes | wall | speedup |
|---|---|---|---|
| 1 | 2 | 19.5 s | 1.00x |
| 2 | 4 | 11.2 s | 1.74x |
| 4 | 8 | 7.4 s | 2.63x |
| 8 | 16 | 6.3 s | 3.07x |
| 11 | 22 | 5.7 s | 3.42x |

3.42x on 11 cores, below linear for reasons that are artefacts of this test
rather than the method: tube and box each spawn `chunks` processes, so
`chunks=11` is 22 processes 2x-oversubscribed on 11 cores; the merge is serial
and included in the timing; and at 20 s total the per-process 0.23 s startup is
a visible fraction. On omega Cen (9.28 CPU-h, 12000 orbits) all three shrink to
noise. Choose `chunks ~ cores/2` since both families run concurrently.

## 8.5 Test 5 -- compiler flags

| variant | chi2 vs baseline | integrate |
|---|---|---|
| baseline (`-ffast-math -O3 -march=native`) | -- | 24.2 s |
| **without `-march=native`** | **identical** | 24.4 s |
| without `-ffast-math` | 0.618% max | 25.8 s |
| neither | 0.618% max | 26.3 s |

**Dropping `-march=native` is free** -- bit-identical chi2, no measurable speed
cost -- and it removes the SIGILL risk on heterogeneous compute nodes. Do it
before the grid. `-ffast-math` is a real choice: it changes results by 0.6%
(comparable to the integration-setting noise in 8.2) and is worth ~7% of
runtime.

## 8.6 Test 6 -- script interactions

`triaxmass`/`triaxmassbin` are gated on `LegacyWeightSolver` and are not
invoked at all under `type: NNLS`. When active they sit once in the tube
branch, and they depend on the potential rather than the orbit library, so a
chunked script must invoke them once, not per chunk. A placement constraint,
not a blocker.

## 8.7 What still needs doing

- Wire chunking into `orblib.py`'s script generation (`chunks` config key,
  N input files, merge step, `triaxmass` invoked once).
- Full omega Cen run as the final gate (9.28 CPU-h).
- Decide the re-baselining question in 8.2. That is a science call: the
  integration-setting sensitivity is real and pre-existing, and chunking
  requires picking one deterministic scheme and living with it.
