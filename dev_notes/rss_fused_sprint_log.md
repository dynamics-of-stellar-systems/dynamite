# RSS/Fused Sprint Log — running learnings (2026-08-21/22)

Companion to `weight_solve_rss_profile.md` (measurements) and
`docs/superpowers/plans/2026-08-21-weight-solve-rss-fused-x.md` (plan).
This file captures what we learned *doing* the work. Update as results land.

## Status board

| item | state |
|---|---|
| Task 1-2: probe + profile driver | done (`ccd6731`, driver `a08249b`) |
| Task 3: baseline measurements | done (`183826f`): peak 626 GiB, ALM resting ~472 GiB median |
| Task 4-7: fused-X path | done + **real-library bitwise proof** (182/182 slab digests identical) |
| Task 8: streamed per-set reads | done, bitwise vs non-streamed on synthetic harness |
| Task 9 step 1: full fused-f64 validation | RUNNING (acceptance: weights bitwise == baseline file; chi2 <=1e-11 rel of 2770835.03357815 / 335126.55535470234) |
| Task 9 steps 2-4: re-profiles x4, concurrency table, docs | queued via `PM_grid/run_task9_queue.sh` |

Commits this leg: ccd6731 b3ab76e 8cc469f 5710209 40fda13 7425ceb 183826f a08249b.

## Key findings

1. **Linux glibc returns huge freed numpy arrays to the OS instantly** -
   reclaim probe: 100% at every size (1-128 GiB), already at `del`, nothing
   left for gc/malloc_trim. The earlier macOS ~0%-reclaim result does NOT
   transfer. Consequence: streamed per-set frees are fully reliable here.
2. **Baseline resting RSS is ~472 GiB median / 626 GiB peak per worker**
   (hists ~165 + A ~124 + X ~124 + ~150 GiB bvls transients). One concurrent
   solve today; fused-X targets X-only residency.
3. **Production shape confirmed**: n_orbs = 45000 (= 2x mirrored tube 30000
   + box 15000); matrix = 371212 x 45000 float64 = ~124 GiB.
4. Fused assembly at production scale: **249 s** to build X (vs ~370 s for
   classic A build) while never allocating A.
5. The ~463 s/iteration ALM pace is bvls-bound and unchanged by this sprint
   by design (chi2-from-resid already landed in the Aug-18 sprint).

## Gotchas that cost us time (do not rediscover)

- **Config paths resolve against CWD**, not the yaml's directory. Every run
  must start in `PM_grid`. The queue script does `cd "$PM"` internally;
  drivers' docstrings now say so.
- **Shell precedence**: `git add ... && git commit ... && nohup job &`
  short-circuits the whole chain when git fails (wrong cwd) - the "launched"
  pid belongs to the dead git job and the log holds a stale traceback.
  Launch long jobs in their own tool call, verify aliveness after ~30 s.
- `_rss_probe` selfcheck race: write one sample row synchronously in
  `__enter__` or the 'init' phase may never appear in the CSV.
- Reclaim probe must FILL the array (`arr[:] = 1.0`); strided touches leave
  pages unmapped and report nonsense (~0.04G for a 64 GiB alloc).
- Binding NNLS staticmethods onto a fake object via `types.MethodType`
  prepends the instance as an extra argument - attach plain functions for
  staticmethods, MethodType only for real methods.
- Stubbed readers must return FRESH copies per call when downstream code
  mutates in place (`scale_x_values` scales x edges in place; a shared stub
  object compounds scaling across calls).
- float32 tolerances are dtype-specific: chi2 identity ~1e-11 rel at f64 but
  ~1e-7 at f32; scaled KKT ratio only ~1e-3 at f32 (cancellation-heavy).
  These are diagnostics, not science outputs; weights remain bitwise.
- `ru_maxrss`/VmHWM is monotonic - phase "peaks" printed from it can be set
  by an earlier transient. Use the CSV's VmRSS column for steady-state
  statements (that is how resting=472 was distinguished from spike=626).

## Open questions for the morning analysis

- Do bvls warm-start transients scale with matrix dtype/size? (drives the
  f32 projections)
- Does streaming measurably lower the ASSEMBLY-window peak below the solve
  window at f64? (if not, its value is f32-only)
- ~~Exact chi2 digits achieved end-to-end vs recorded baseline~~ ANSWERED,
  see gate revision below.

## GATE REVISION (2026-08-22, after full fused validation)

The planned acceptance "weights bitwise == recorded baseline file" was
**mis-calibrated**: adelie's parallel BVLS is not bitwise-reproducible
across processes (thread-scheduling reduction order), and the original
baseline additionally predates `6fe911c`, whose chi2-from-resid change can
flip best-iterate selection among near-tied multiplier updates.

Measured end-to-end (fused-f64 full run vs recorded original, both 100 ALM
iterations, identical config):

| metric | result |
|---|---|
| peak RSS | **503.0 vs 626.1 GiB (-123 GiB = A eliminated)** |
| chi2_tot | 2770837.5186 vs 2770835.0336 -> rel diff 9.0e-07 |
| chi2_kin | 335126.6108 vs 335126.5554 -> rel diff 1.7e-07 |
| weights | max|dw| 2.8e-06; support differs by 13 of 45000 orbits; sum(w) agrees to 2.4e-9 |

Verdict: statistically identical solutions of an identically-posed problem.
What stays bitwise-proven: X/col_norm/y inputs (real-library digest check)
and construction-path equality; what was never well-posed: bitwise output
equality across processes/versions of a threaded solver.

Revised acceptance (met): chi2 within <=1e-06 relative; support overlap
>99.8%; sum(w) gap unchanged (~4e-8, the ALM constraint tolerance).
