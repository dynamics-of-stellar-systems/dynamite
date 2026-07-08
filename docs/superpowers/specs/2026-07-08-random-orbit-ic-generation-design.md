# Random (Agama-backed) Orbit Initial-Condition Generation — Design

## Motivation

Dynamite currently seeds Schwarzschild orbit libraries from a regular grid in
the integrals of motion (E, I2, I3), via the Fortran `orbitstart`/
`orbitstart_bar` binaries. Vasiliev (2012) showed this discrete grid sampling
can produce artifacts in the resulting model. Both `Forstand`
(Vasiliev 2020, built on `Agama`) and the newer `SchwarMAX`
(Zhang et al. 2026, arXiv:2607.00117) instead sample initial *positions* from
a smooth fiducial density and initial *velocities* from a local
Jeans-equation velocity ellipsoid, producing a smoother, more
equilibrium-like seed with only light dithering, rather than heavy grid
dithering.

This project explores whether the same approach is worth adopting in
dynamite, and if so, does it with the smallest possible footprint on
dynamite's existing (mostly Fortran) orbit-library pipeline.

## Goal

Add a new, opt-in orbit initial-condition (IC) generation path that:
- Samples orbit seed positions from a smooth fiducial density and velocities
  from Agama's Jeans-equation solve, instead of a regular (E, I2, I3) grid.
- Requires **zero changes** to the Fortran orbit integrator, orbit
  classification, mirroring, or NNLS weight-solving — every downstream
  consumer of the orbit library is unaffected.
- Is a standalone feature, independent of the (unfinished, possibly
  abandoned) `jax-orblib` branch's JAX orbit-integrator work.

## Scope

**In scope:**
- Non-rotating (static/triaxial) systems only — matches dynamite's
  primary current use case (e.g. NGC5139/ω Cen). Bar-disk/rotating
  (Jacobi energy) systems are not addressed.
- Keeping dynamite's existing tube/box orbit-family split for
  integration purposes (see "Why keep the tube/box split" below).
- A new Python module plus a config flag; validated on a small existing
  `dev_tests` config, not a full ω Cen production run.

**Out of scope:**
- Bar-disk/rotating systems.
- Any change to `jax-orblib` or reliance on it landing.
- Reproducing the Fortran `noreg` field's original intent — confirmed
  dead code (see below), not carried forward at all.
- Reproducing the Fortran integrator, classification, mirroring, LOSVD/
  qgrid accumulation, or weight-solving in any new backend (Agama or
  otherwise) — approach considered and explicitly rejected as
  substantially more work (see "Alternatives considered").

## Why keep the tube/box split

Investigated whether dynamite's a priori tube/box orbit-family split
(separate `orbitstart`/`orbitstart_bar` seeding, separate `orblib`/
`orblibbox` Fortran integration, separate mirroring rules) is something a
continuous/random sampler would need to reproduce or restructure.

Confirmed it does not: `integrator_find_orbtype`
(`legacy_fortran/orblib_f_new_mirror.f90:757-828`) classifies every
integrated orbit's actual dynamical type (X/Y/Z-tube, box, chaotic) purely
from the sign-consistency of Lx/Ly/Lz along its trajectory — computed
identically and unconditionally regardless of whether the orbit came from
`begin.dat` or `beginbox.dat`. The IC-time tube/box split is a seeding/
mirroring convenience only, not the actual physical classification, which
is redone from scratch post-integration.

This means we can keep feeding the existing two Fortran integrator
binaries unchanged: seed a "tube-family" subset with full 3D velocities
from the Jeans ellipsoid, and a "box-family" subset with velocities forced
to `v_y=0` (matching the existing a priori box-seeding convention), and
write them to `begin.dat`/`beginbox.dat` respectively.

## Why the file interface is a safe seam

Confirmed the Fortran integrator (`integrator_setup()`/
`integrator_setup_bar()`, `orblib_f_new_mirror.f90:206-252`) reads orbits by
sequential file position, grouping every `dithering³` consecutive lines
into one orbit bundle, using only the aggregate counts (`nEner*nI2*nI3`,
`orbit_dithering` from config) for that grouping — not the per-line
`iE`/`iI2`/`iI3` values. So a replacement IC generator only needs to:
- write the same fixed-width Fortran format (`(3I5, 9ES30.10, I4)`,
  documented in `dev_notes` from the `jax_orbit_integration` plan and
  verified against `orbitstart_f.f90`),
- group lines sequentially into blocks of `dithering³`,
- and set header integers (`nE nI2 nI3`) whose product equals the total
  line count.

The per-line grid-index fields are along for the ride, not load-bearing.

## `rcirc`/`vcirc`/`tcirc`

These come from simple circular-orbit evaluations against the potential
(`orbitstart_f.f90:83-86`, via `ip_accel`/`ip_potent`) — not tied to the
grid. Agama's `circularVelocity(R)` gives this directly against an Agama
`Potential` object, so no custom implementation needed.

## `noreg` is dead code — dropped entirely

`noreg` ("not regularizable") was originally a numerical-stability flag:
`find_unregorbits` (`orbitstart_f.f90:163-181`) marks orbits at grid
boundaries between regular/irregular energy shells or at the
tube/box family boundary, so the integrator could skip orbit
regularization (a coordinate-transform trick used near degenerate/
central orbits) for those borderline cases.

Checked whether this field is read anywhere downstream: it is not.
Neither `orblib_f_new_mirror.f90` (the active integrator) nor the older
`orblib_f.f90`, nor the Python `orblib.py`, reads the `noreg` column back
from the IC file. It is written but never consumed — dead code, likely a
leftover from an earlier integrator version. The new IC generator will
simply not emit meaningful values for this field (write a constant `0`)
and this is not treated as a simplification/gap, since nothing downstream
would differ either way.

## Architecture

A new module, `dynamite/random_ic.py`, is invoked at the exact point
`orblib.py` currently shells out to the `orbitstart`/`orbitstart_bar`
Fortran binary. A new config option, `orblib_settings.ic_generator: grid |
random` (default `grid`), selects between the existing path and this one.
When `random` is selected, `dynamite/random_ic.py` runs instead of the
`orbitstart` subprocess call, and writes `begin.dat`/`beginbox.dat` in the
same format the Fortran binary would have. Everything downstream is
unmodified and unaware of which path produced its input files.

## Data flow

1. Build an Agama `Potential` object from the trial model's MGE + dark
   halo + BH parameters (the config already carries all of these — this
   is a translation layer, not new physics).
2. Sample `N_bundles` positions from a fiducial tracer density using
   Agama's sampling utilities.
3. For each position, get the local velocity ellipsoid from Agama's
   Jeans-equation solve. Draw one full-3D velocity for the "tube-family"
   subset; for a separate "box-family" subset, force `v_y=0` (matching
   the existing a priori box-seeding convention).
4. Dither: generate `dithering³` nearby realizations per seed, reusing
   the existing `dithering` config value, reinterpreted as
   "realizations per seed" rather than "grid-cell sub-sampling".
5. Compute `rcirc`/`vcirc`/`tcirc` per seed via Agama's
   `circularVelocity`.
6. Write `begin.dat` (tube) and `beginbox.dat` (box) sequentially,
   grouped in blocks of `dithering³`, with header integers whose product
   equals the total count.

## Key components

- `dynamite/random_ic.py` — potential translation, sampling, dithering,
  file writing (new).
- `dynamite/config_reader.py` — add `ic_generator` setting (default
  `grid`, validated against `{grid, random}`).
- `dynamite/orblib.py` — at the current `orbitstart` subprocess call
  site, branch on `ic_generator` to call `random_ic` instead.
- A new dev test config — a copy of `dev_tests/user_test_config_ml.yaml`
  with `ic_generator: random`, for direct A/B comparison against the
  existing grid-based run on the same tiny system (nE=6, nI2=5, nI3=4,
  dithering=1 → 120 orbits).

## Testing / validation

Run the small dev_tests config twice (grid vs. random), same potential,
and compare:
- orbit-library-level sanity: energy conservation per orbit, a reasonable
  spread of orbit classifications from `integrator_find_orbtype` (not
  everything collapsing to "chaotic"),
- final NNLS fit quality (chi2) on whatever mock/reference data that
  config already uses.

This is a proof-of-concept check on a tiny system, not a rigorous
benchmark against ω Cen.

## Alternatives considered

**Full Agama-based orbit integration** (not just IC generation) —
rejected. Agama's orbit integrator itself is solid and would be free, but
this would discard everything dynamite's Fortran currently does
downstream of integration: LOSVD histogram accumulation, 3D intrinsic-
mass (`qgrid`) accumulation, sky-plane projection + PSF convolution,
orbit classification and tube/box mirroring, and dithering/bundle
bookkeeping for NNLS. All of that would need rebuilding against Agama's
orbit output format — a bigger scope than the still-unfinished
`jax-orblib` plan, which only tackles the integrator/accumulator and
still leaves classification/mirroring in Fortran.

**Grid-jitter without Agama** — a cheaper first increment considered:
perturb the (E, I2, I3) grid values with quasi-random offsets in Python,
feeding the existing unmodified `orbitstart_f.f90`. Would test whether
avoiding exact grid points alone helps, without building the Agama/Jeans
machinery. Not pursued as the primary path since the goal is the real
Forstand-style equilibrium seeding, but noted as a possible cheap
sanity-check step if useful later.

## Open questions / risks for the implementation plan

- Agama installability in the target environments (this machine, the
  `gc-diag` conda env, the cluster `ENV`) is unverified — first
  implementation task should confirm this before any other work.
- Exact on-disk fixed-width format must be re-verified byte-for-byte
  against a real `orbitstart`-produced `begin.dat`, not just the
  `jax_orbit_integration` plan's documentation of it.
- Box-family fraction (what proportion of seeds get `v_y=0` forced) is
  an open tuning parameter, not yet specified.
