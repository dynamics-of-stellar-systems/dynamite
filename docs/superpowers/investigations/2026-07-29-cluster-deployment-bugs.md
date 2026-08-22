# Cluster Deployment: Affinity Bug, Redundant Builds, and Chunk Merge

Date: 2026-07-29
Branch: `feature/orblib-performance`
Config: NGC5139 (ωCen), nE=40 nI2=20 nI3=15 dithering=3, 12000 starting points
Node: Node-12 (192 logical cores, non-exclusive interactive session)

## Problem 1: CPU Affinity Pinned to Core 0

**Symptom:** All orblib and weight-solver processes showed ~1% CPU despite correct
config sizing (total_cores=90, orblib_chunks=auto→5). `taskset -pc <pid>` showed
affinity list: `0` for every Fortran and Python worker process.

**Root cause:** `import adelie.solver` (the Rust/rayon-backed NNLS solver) triggers
rayon's thread-pool initialization, which mis-detects the core count on this
192-core machine and narrows the *importing process's own* CPU affinity to a
single core as a side effect. Every descendant (forked pool workers, subprocess
chunks) inherits that single-core mask.

**Diagnosis path:**
1. Ruled out OMP/BLAS threading vars (OMP_NUM_THREADS, OPENBLAS_NUM_THREADS,
   MKL_NUM_THREADS) — verified they were set to 1, not the cause.
2. Ruled out OMP_PROC_BIND — test harness showed both TRUE and FALSE gave full
   0-191 affinity in a plain-numpy subprocess test.
3. Ruled out pathos multiprocessing pool — affinity_test2.py showed pool workers
   and their subprocess children retain full 0-191.
4. Ruled out sustained-load cgroup throttling — affinity_test3.py ran 40 CPU-bound
   workers for 2 min, all stayed at 192 CPUs.
5. Ruled out nested subprocess pattern — affinity_test4.py matched the real
   chunked bash-script-with-background-children structure, still full affinity.
6. **Walked the actual process tree** with `taskset -pc` at every level:
   - tmux → bash: 0-191
   - python oCen_grid.py: **0** ← narrows here
   - bash script → orblib_new_mirror: 0 (inherited)
7. **Bisected imports** using `os.sched_getaffinity` after each submodule import
   in `dynamite/__init__.py` → collapsed at `import dynamite.weight_solvers`.
8. **Isolated to `import adelie.solver` alone** — confirmed `sched_getaffinity`
   drops from 192→1 on that single import.

**Fix** (`dynamite/weight_solvers.py`, commit 7459b21 → amended into d6b1712):
```python
if hasattr(os, 'sched_setaffinity'):
    os.sched_setaffinity(0, range(os.cpu_count()))
```
Placed immediately after `import adelie.solver`, before any forking.

**Verification:** `os.sched_setaffinity(0, range(os.cpu_count()))` restores 192.
Confirmed end-to-end: grid ran at 100% CPU per process after reinstall.

---

## Problem 2: Redundant Orbit Library Integration (4× per library)

**Symptom:** Each of 9 orbit libraries was built 4 times concurrently (log showed
4 "Integrating orbit library for ..." per directory), producing 271 orblib
processes against 90 total_cores budget.

**Root cause:** `get_orblib()` has a racy check-then-build:
```python
if not os.path.isfile(self.mod_dir + 'datfil/tube_box_done'):
    # ... build
```
Multiple concurrent models sharing the same orbit library all see
`tube_box_done` missing and all start building. The existing code only cleaned up
*after* redundant builds happened ("already merged by a concurrent model,
discarding this process's chunks"), it didn't prevent them.

**Fix** (`dynamite/orblib.py`, commit d6b1712):
- `claim_orblib_build()` — atomic lock file via `os.open(path,
  os.O_CREAT|os.O_EXCL)` with stale-lock recovery (checks if the locking process
  still exists via `os.kill(pid, 0)`).
- `wait_for_orblib_build()` — polls for `tube_box_done` with timeout,
  falls back to building if lock goes stale.
- `release_orblib_build_claim()` — removes lock file in `finally` block after
  build completes (whether success or failure).

**Verification:** `grep "claim" dynamite.log` should show exactly 9 lines (one
per library) — zero claim lines in a run means the old code path is still
active, likely from stale `tube_box_done` files surviving across runs.

---

## Problem 3: Chunk Merge Writes Wrong Orbit Count (11988 ≠ 12000)

> **RESOLVED 2026-07-30 — the fix recorded below was wrong and made things
> worse.** Patching the header hid the real fault (the chunks only covered
> 1/dithering³ of the library) and converted a correct error message into the
> unintelligible crash recorded as Problem 4. See "Resolution" at the end of
> this document. The header is now checked, not overwritten, and the body's
> orbit count is checked too.

**Symptom:** "Number of orbits in ... is 11988, but expected 12000." Every
library, both families, every iteration.

**Root cause:** The chunked merge in `get_orbit_library_chunked` computed
`n_orbits` as `n_orbit_starting_points() * dithering**3`, which for this config
is `(40*20*15 // 27) * 27 = 444 * 27 = 11988`. This value was written into the
qgrid file header (byte 4). But `read_orbit_base` (the reader) expects
`nE * nI2 * nI3 = 40 * 20 * 15 = 12000`, matching what the unchunked Fortran
output always contains.

**Fix** (amended into d6b1712):
```python
# Before (wrong — post-dithering count):
n_orbits = self.n_orbit_starting_points() * self.settings['dithering'] ** 3
# After (correct — matches reader):
n_orbits = self.settings['nE'] * self.settings['nI2'] * self.settings['nI3']
```

Note: `n_orbit_starting_points()` is still correct for `orbit_chunk_bounds()`
(chunking needs starting-point ranges). Only the header-patch value was wrong.

**Verification on cluster:** All 9 libraries report 12000 orbits.

---

## Problem 4: Read Corruption — All Paths (chunked AND unchunked)

> **RESOLVED 2026-07-30.** The hypothesis below (a newer Fortran binary writing
> a different file format) was wrong, and so were two of the observations it
> rested on. See "Resolution" at the end of this document.

**Symptom:** Read crashes with "Size obtained (1) is not a multiple of the dtypes
given (4)" in `_read_individual_orbit` (orblib.py:1269). The scipy FortranFile
encounters a record of 1 byte where it expects multiple int32. This happens on
EVERY model, every iteration, persistently.

The qgrid file header says 12000 orbits, but walking the raw records finds only
1332 before the trailing 1-byte record — the rest of the file is not valid
Fortran records. Similarly, the losvd file header record starts with (1012, 66,
...) instead of (12000, 1012, ...), and walking it finds 450k records — the
content is there but incorrectly structured.

**Crucial finding: The crash reproduces with `orblib_chunks: 1` (unchunked
path).** This rules out the chunk merge, the Fortran binary's range handling,
and the n_orbits header patch as causes. The issue is in the unchunked Fortran
output itself — the `orblib_new_mirror` binary on this cluster produces files
that `read_orbit_base` cannot parse.

**Hypothesis:** The installed Fortran binary is a newer version than the Python
reader expects. The binary may write a different file format (e.g., different
record layout, additional leading record) that the old version of
`_read_individual_orbit` doesn't handle. This would explain why the record walk
finds data beyond the header but with wrong record sizes — the file's actual
structure is different from what the reader assumes.

**Next steps:**
1. Check the binary's timestamp vs the repo's source:
   ```
   ls -la /nexus/.../legacy_fortran/orblib_new_mirror
   ```
2. Compare a small test file's record layout (e.g., dump the first ~1000 bytes
   as hex) against what `_read_individual_orbit` expects.
3. If the format differs, either rebuild the binary from the repo's Fortran
   source, or update the Python reader to match the installed binary's format.

---

## Cluster Environment Notes

- Python: conda env at `/nexus/posix0/MIA-astro-env/nneum/pesmith/ENV/`
- Dynamite installed in site-packages (NOT repo checkout — repo at ~/research/dynamite)
- Fortran binaries in:
  `/nexus/posix0/MIA-astro-env/nneum/pesmith/ENV/lib/python3.12/site-packages/legacy_fortran/`
- To reinstall after pulling: `pip install -e /path/to/repo --force-reinstall --no-deps`
- Node-12 is interactive (shared), not batch-scheduled — jsauter's ft_fit.py
  processes occupied ~19 cores throughout.

## Key Commands for Diagnosis

```bash
# Affinity check
taskset -pc <pid>
python3 -c "import os; print(len(os.sched_getaffinity(0)))"

# Read orbit count from bz2 header
python3 -c "
import bz2, struct
raw = bz2.open('NGC5139_output/models/*/datfil/orblib_qgrid.dat.bz2').read()
n = struct.unpack('<i', raw[4:8])[0]
print(n, 'orbits')
"

# Walk orbit records in a qgrid file
python3 -c "
import bz2, struct
raw = bz2.open(f).read()
pos = 0; orb = 0
while pos < len(raw) - 8:
    lead = struct.unpack('<i', raw[pos:pos+4])[0]
    if lead == 1: break
    trail = struct.unpack('<i', raw[pos+4+lead:pos+8+lead])[0]
    if lead != trail: break
    orb += 1; pos += 8 + lead
print(orb, 'records')
"
```

---

# Resolution (2026-07-30)

Problems 3 and 4 were one bug, and the "fix" for Problem 3 caused Problem 4.
The Fortran binary was never at fault.

## Root cause

`n_orbit_starting_points()` in `dynamite/orblib.py` returned

```python
(nE * nI2 * nI3) // dithering**3        # 444 for omega Cen
```

but the legacy Fortran has already applied that factor by the time it counts
starting points:

- `iniparam_f.f90:159-161` multiplies the grid read from `parameters_pot.in`
  by `orbit_dithering` (`Nener = Nener*orbit_dithering`, same for nI2/nI3).
- `orbitstart` writes *that* dithered grid into `begin.dat`'s header.
- `integrator_setup_write` divides by `dithering**3` again to get the bundle
  count.

The two cancel: the number of starting points — what the Fortran's
`[starting orbit]` and `[orbits to integrate]` inputs count, and therefore what
a chunk range refers to — is `nE * nI2 * nI3` = **12000**, not 444.

Verified empirically against `orbitstart` at both dithering values:

| config | `begin.dat` header | entries | starting points |
|---|---|---|---|
| 6,5,4 d=1 | `6 5 4` | 120 | 120 |
| 6,5,4 d=2 | `12 10 8` | 960 | 960/8 = 120 |
| 40,20,15 d=3 (ωCen) | `120 60 45` | 324000 | 324000/27 = 12000 |

Independently confirmed by the Fortran's own bound: `integrator_number` is
capped at `(nEner*nI2*nI3/integrator_dithering**3)` with `stop " Too many
orbits in total"`, so chunk ranges summing to 12000 are exactly the maximum
allowed. Ranges are 1-based inclusive, matching `orbit_chunk_bounds`.

## Why the chain of wrong conclusions

1. **The chunks covered 444 of 12000 starting points**, so the merged qgrid held
   444 orbits = 1332 records (3 per orbit, from `_read_individual_orbit`), then
   the closing 1-byte record. That is exactly the record count in Problem 4.
2. **The natural merged header value was 444×27 = 11988** — the Problem 3
   symptom. Patching it to 12000 silenced the one check that was correctly
   reporting the shortfall, so the reader looped 12000 times over a 444-orbit
   body, ran off the end, and hit the closing record: *"Size obtained (1) is not
   a multiple of the dtypes given (4)."*
3. **The losvd header `(1012, 66, ...)` was never corrupt.** It holds
   `(nconstr, nvcube, dvcube)` — see `histogram_setup_write` and the reader's own
   comment at `orblib.py:1438-1442`. It has never contained an orbit count.
4. **`orblib_chunks: 1` did not really reproduce it.** That path cannot produce
   a short library. It was almost certainly reading a truncated `.bz2`: the
   parallel unchunked script rewrote the archive in place
   (`rm -f f.bz2 && bzip2 -k f`), so a concurrent model decompressed a partly
   written file — which reads as a good header followed by garbage records, the
   same signature as a short library.

Only the chunked path broke because the unchunked path passes
`number_orbits: -1`, which the Fortran expands to the full 12000 itself.

## Why it was never caught

`dev_tests/user_test_config_ml.yaml`, the configuration the chunking gate ran
with, uses `dithering: 1`, where `// dithering**3` is the identity. ωCen at
`dithering: 3` was the first configuration that could expose it. The gate also
compared only the four binary library files (so a missing `orbclass` file passed
unnoticed) and reused a pre-existing model tree's `begin.dat`, which at one
point meant it ran against initial conditions for a different grid entirely.

## Other bugs found in the same code

All fixed in the same series of commits:

- **`orblib.dat_orbclass.out` was destroyed by chunking** — chunks deleted,
  never merged; read by `read_orbit_property_file`, `analysis.py` and
  `plotter.py` for λ_z, and counted by `model.py` towards a complete library.
- **`tube_done`/`box_done` never created after a chunked build**, so
  `tube_box_done` never appeared and every model reintegrated the library.
- **Non-atomic `.bz2` publish** in the parallel unchunked script (item 4 above).
- **`triaxmass` ran against a non-existent merged library** under
  `LegacyWeightSolver` + chunking; chunking is now refused for it, and for
  libraries restricted via `starting_orbit`/`number_orbits`.
- **2× core oversubscription** in `orblib_chunks: auto` when
  `orblibs_in_parallel: False` — the chunked script always runs both families.
- **A failed read left the pool worker in the wrong directory**, so unrelated
  models afterwards resolved `datfil/...` against the previous model's path.
- **A failed read leaked the decompressed library** (~2 GB per family at ωCen
  scale) on every failed model of every iteration.
- **Stale `begin.dat` silently reused** when the orbit grid changed, since model
  directories are `orblib_<iteration>_<row>` and not keyed on the grid. This
  produces the Fortran's opaque `STOP  Not so many orbits`, whose bound comes
  from `begin.dat` rather than the configuration.
- **`wait_for_orblib_build` rebuilt a completed library** when the owner's
  done-file and lock-release fell inside one poll interval, and recursed once
  per hour when the build legitimately took longer than the timeout.
- **Unintelligible failures**: `bunzip2` exit codes were ignored, and a
  truncated losvd raised a bare `IndexError` from inside the record walk.

## Verification

- `dev_tests/test_orbit_chunk_bounds.py` — chunk ranges cover `nE*nI2*nI3` and
  are consecutive; fails on the old arithmetic.
- `dev_tests/test_orblib_chunk_merge.py` — synthetic files at the exact record
  layout; byte-identity across five chunk splits, the vectorised reader walking
  the merged losvd, and rejection of both a wrong header and a short body.
- `dev_tests/test_orblib_failure_hygiene.py` — cwd and temp-file hygiene on a
  failed read, named error on a truncated losvd, stale-ICS detection for all
  four grid settings, and `can_chunk_orbits` refusing the unsafe cases.
- `dev_tests/test_orblib_chunking.py` end to end at **`dithering: 2`** (chunks
  1/2/3/5/7) and `dithering: 1` (chunks 1/4/8): byte-identical libraries,
  including both `orbclass` files, and identical chi2/kinchi2. The reverted
  arithmetic reproduces the cluster error locally.

## Before rerunning on the cluster

The libraries already on disk are genuinely short, and nothing will rebuild them
while they look complete:

```bash
rm -f NGC5139_output/models/*/datfil/{orblib,orblibbox}_{qgrid,losvd_hist}.dat.bz2 \
      NGC5139_output/models/*/datfil/{orblib,orblibbox}.dat_orbclass.out \
      NGC5139_output/models/*/datfil/{tube,box,tube_box}_done \
      NGC5139_output/models/*/datfil/building.lock \
      NGC5139_output/models/*/datfil/*_qgrid.dat.tmp
```

## Known remaining limitations

- `merge_files` loads every chunk into memory: ~2.5 GB peak per family at ωCen
  scale, one family at a time. Measured, and small against the ~100 GB per-model
  read peak, so left alone; streaming the copy is the fix if it ever matters.
- The losvd stream has no independent orbit count of its own — `orblib_chunks`
  does not know the aperture count — so it relies on the qgrid check, which
  covers the same chunk ranges.
- Per-chunk `infil/*.in` files and `cmd_tube_box_orbs_chunked_p<pid>` scripts are
  never cleaned up. Harmless clutter.
