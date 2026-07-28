# Parallelising orbit integration by splitting the orbit range

Status: **planning / feasibility**. Created 2026-07-28.

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
