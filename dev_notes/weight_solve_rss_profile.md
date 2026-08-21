# Weight-Solve RSS Profile — Baseline Measurements on the Cluster Node

Date: 2026-08-21 · Machine: 192-core / 1416 GB KVM node, NFS scratch
Branch: `feature/orblib-performance` @ fused-sprint worktree (Tasks 1-8 of
`docs/superpowers/plans/2026-08-21-weight-solve-rss-fused-x.md`)
Run: `_real_solve_rss_profile.py` on the real xeast orblib (nE×nI2×nI3 =
30×25×20, tube+box mirrored → **n_orbs = 45000**; matrix 371212 × 45000
float64 ≈ **124 GiB** per full matrix), adelie solver capped at
`--alm-iters 30`, OMP_NUM_THREADS=24. Timeline CSVs:
`PM_grid/rss_base_f64_cap30.csv` (+ `_summary.json`).

## Headline answers (spec Phase-1 questions)

1. **Histograms weigh far more than extrapolated**: orblib read peaks at
   **164.4 GiB** (the old ~34.7 GB figure was laptop-scale). All kinematic
   sets' combined histograms are retained through assembly and solve today.
2. **Linux glibc reclaims freed huge allocations immediately** (reclaim
   probe below: 100% at every size, already at `del`, before gc/malloc_trim).
   The macOS ~0%-reclaim result does not apply here → streamed frees are
   fully reliable.
3. **Resting RSS during the ALM loop is ~472 GiB median** (floor 377 GiB,
   repeating transients to ~626 GiB). Global peak **626.1 GiB**. This is the
   number that caps concurrency at one solve today.
4. Combine/mirror transient is folded into the read phase (peak 164 GiB);
   no separate spike beyond it was observed.

## Reclaim probe (`_rss_probe.py --sizes 1 8 64 128`)

```
   GiB   rss@alloc   after_del   after_gc  after_trim reclaimed%
   1.0       1.04G       0.04G       0.04G       0.04G    100.0%
   8.0       8.04G       0.04G       0.04G       0.04G    100.0%
  64.0      64.04G       0.04G       0.04G       0.04G    100.0%
 128.0     128.04G       0.04G       0.04G       0.04G    100.0%
```
(after_del = after_gc = after_trim: numpy's mmap-backed buffers are unmapped
at `del` on this glibc; `malloc_trim` has nothing left to do.)

## Phase-by-phase timeline (VmRSS, GiB)

| phase | samples | min | median | max |
|---|---|---|---|---|
| orblib_read | 1109 | 0.4 | 15.8 | 164.5 |
| matrix_build | 156 | 129.0 | 208.0 | 252.2 |
| alm_solve:enter→iter_1 | 234 | 252.2 | 368.9 | 501.1 |
| alm iterations (per-10 blocks) | ~900 | 376.7 | 472–503 | 625.7 |
| kkt_violation | 18 | 377.2 | 421.2 | 501.7 |
| post-solve | 1 | — | — | 252.8 |

Durations (summary JSON): wall 10195 s ≈ 2.83 h total; read 1428 s;
matrix build 369 s; ALM ≈ 463 s/iteration at 24 threads.

Interpretation:
- Steady solve residency ≈ hists (~165) + A (~124) + X (~124) plus solver
  temporaries → the ~377 GiB floor; each bvls call transiently adds ~150 GiB
  (warm-start state + internal copies), spiking HWM to 626.
- The X-build moment itself (A+X+hists all alive) coincides with the first
  iteration's spike.
- Post-solve the process drops to ~253 GiB (hists + weights), confirming
  that freeing works cleanly once references die.

## Numerics recorded by this run

chi2_tot = 2871744.4236285393, chi2_kin = 335125.8343159734 at 30 ALM
iterations (capped run - NOT converged to the full baseline's
2770835.03357815; chi2_kin is already close because the kinematic terms
converge early). Driver-saved weights are bitwise equal to the
`orbit_weights.ecsv` written by solve() (verified).

## Per-worker budget & streaming go/no-go

| mode | expected steady peak | fits in 1416 GB (0.85 safety) |
|---|---|---|
| baseline classic (today) | ~500-626 GiB | 1 worker |
| fused-X float64 (X only) | ~124 GiB + one-set hists | 8 workers |
| fused-X float32 | ~62 GiB + one-set hists | up to 16 workers |

Streaming saves the difference between *all* sets' histograms (~165 GiB)
and the largest single set during the assembly window only; with fused-X the
assembly window is no longer the global peak, so streaming's value is
bounded but still material for float32 (where one PM set could rival X).
Proceed with streamed reads as planned; measure both in the re-profiles.
