# Implementing chunked orbit integration in DYNAMITE

Status: **planned**. Created 2026-07-28.
Feasibility study and measurements: `2026-07-28-orbit-integration-parallelisation.md`.

The Fortran side is already done and committed on this branch: an orbit library
can be built by N processes over disjoint orbit ranges, and the result is
bit-identical to a single-process run. What remains is exposing that through
`orblib.py` and the configuration file, documenting it, and validating it at
full scale.

---

## 1. Goal and acceptance criteria

**Goal.** Let a user split each orbit family across N processes with one
configuration setting, and have the result be indistinguishable from today's
single-process library.

**Acceptance.** Two full-scale omega Cen models, identical in every respect
except the new setting, must produce:

1. **byte-identical** `orblib_qgrid.dat.bz2`, `orblib_losvd_hist.dat.bz2` and
   the box equivalents, after merging;
2. **identical** `chi2` and `kinchi2` from the weight solver;
3. a wall-clock reduction consistent with the chunk count.

(1) is the strong claim and the one worth designing the whole exercise around.
It already holds at NGC6278 scale for 4 and 8 chunks; section 5 is about
demonstrating it where it matters.

## 2. Configuration interface

### 2.1 Where the setting belongs

In `multiprocessing_settings`, beside `ncpus` and `orblibs_in_parallel`::

    multiprocessing_settings:
        ncpus: 192
        ncpus_weights: 40
        orblibs_in_parallel: True
        orblib_chunks: 4          # NEW: split each orbit family across N processes

This is a resource decision, not a physical one: it cannot change results, so
it does not belong in `orblib_settings` next to `nE`/`nI2`/`nI3`, which define
the library itself.

Note this placement does **not** avoid the config-diff warning.
`Model.validate_config_file` diffs the whole file line by line, so changing
`orblib_chunks` will raise "ACTION REQUIRED: the current config file differs
from the config file backup" against every existing model, even though the
setting provably cannot alter them. That is pre-existing behaviour and out of
scope here, but it must be called out in the documentation, and it shapes the
validation protocol in section 5 (use separate output trees).

### 2.2 Semantics

- Type: integer >= 1. **Default 1**, which reproduces today's behaviour exactly.
- Applies to each orbit family independently. Processes per model are therefore

      orblib_chunks x (2 if orblibs_in_parallel else 1)

  which extends the process-budget rule already documented for
  `multiprocessing_settings`:

      ncpus x orblib_chunks x families x OMP_NUM_THREADS  <=  N_CPU

- Orbits are divided as evenly as possible, with the remainder distributed one
  per chunk rather than all landing on the last one.
- `orblib_chunks: 1` must take the existing code path unchanged, not a
  one-chunk special case of the new one. This keeps the default risk-free and
  makes the A/B in section 5 meaningful.

### 2.3 Validation, in `config_reader.py`

Alongside the existing `orblibs_in_parallel` handling (~line 537):

- default to 1 when absent; reject non-integer or < 1 with a clear message;
- clamp to the orbit count `nE*nI2*nI3 / dithering**3` with a warning, since
  more chunks than orbits is a user error, not a crash;
- log the resulting process count per model, so an over-subscribed run is
  visible in the log rather than only in the wall clock;
- warn (do not fail) when `ncpus x orblib_chunks x families` exceeds the
  detected CPU count.

An `orblib_chunks: 'auto'` value, deriving the count from spare cores, is
deliberately **not** in scope. It needs to know how many models run
concurrently, which lives in the model iterator, and a fixed integer gets the
entire benefit for a cluster run. Revisit once the search driver is settled.

## 3. Code changes

All in `dynamite/orblib.py` unless noted. Roughly 100 lines.

| # | Change | Notes |
|---|---|---|
| 1 | `write_orblib_dot_in(box)` -> `(box, start=None, number=None, tag=None)` | `starting_orbit`/`number_orbits` are already config keys threaded into this function, and every output filename derives from the single `o_file` variable, so tagging is one line |
| 2 | `write_executable_for_integrate_orbits_chunked()` | Model on `..._par()`, which already backgrounds jobs and `wait`s. Emits `2 x chunks` invocations |
| 3 | `get_orbit_library_chunked()` | Copy of `get_orbit_library_par()` plus a merge call once the script returns |
| 4 | Move the merger into the package | `dev_tests/orblib_merge_chunks.py` -> `dynamite/orblib_chunks.py` |
| 5 | Dispatch in `get_orblib()` (~line 138) | Extend the existing `orblibs_in_parallel` branch |
| 6 | `triaxmass`/`triaxmassbin` invoked once | After the chunk `wait`, not per chunk. Gated on `LegacyWeightSolver`, so unused under `NNLS`, but must be correct |
| 7 | `config_reader.py` validation | Section 2.3 |

Two details that will bite if forgotten, both learned the hard way:

- **Each chunk writes a trailing 1-byte Fortran record.** A merged file must
  carry exactly one of them, not one per chunk. The merger handles this by
  locating the final record from the trailing length marker.
- **The `.tmp` resume file is keyed on the output filename.** Distinct chunk
  tags therefore keep chunks from resuming each other's work, but the tags must
  actually differ per chunk *and* per family.

## 4. Documentation

- **Docstrings** on every new method, numpydoc style to match the module:
  parameters, returns, raises. `get_orbit_library_chunked` should state that
  the merged output is bit-identical to the unchunked path, since that is the
  property a reader will want to rely on.
- **`docs/getting_started/configuration.rst`**, `multiprocessing_settings`
  section: document `orblib_chunks` next to `orblibs_in_parallel`, and extend
  the existing `OMP_NUM_THREADS` recommendation to include the chunk factor.
  That section already reasons about a process budget, so this is an extension
  of an existing explanation rather than a new one.
- **A short note on when it helps.** Chunking reduces wall clock, not CPU. On a
  saturated node it buys nothing; its value is latency, when a search iteration
  has fewer models than the node has cores. Users will otherwise enable it and
  be puzzled that throughput did not improve.
- **The config-diff warning** (section 2.1), so the first person to change the
  setting on an existing tree is not alarmed.

## 5. Full-scale A/B validation

Two omega Cen runs, identical except `orblib_chunks`.

**Setup.** Separate output trees via `io_settings`, so each keeps its own
config backup and the diff warning in 2.1 does not confuse the comparison.
Same `random_seed`, same `accuracy`, same everything else. Run A with
`orblib_chunks: 1`, run B with `orblib_chunks: 8` (16 processes with
`orblibs_in_parallel: True`).

**Compare.**

| Check | Expectation |
|---|---|
| `orblib*_qgrid.dat.bz2`, `orblib*_losvd_hist.dat.bz2` | byte-identical (compare uncompressed; bzip2 output should also match, but the uncompressed files are the claim) |
| `chi2`, `kinchi2` | identical |
| orbit count in the merged qgrid header | 12000 per family |
| wall clock | B substantially below A's 5.87 h |

**A failure here means one of two things**, and they are worth distinguishing:
either the merge is wrong (a file-format problem, which the byte comparison
localises immediately), or another piece of order-dependent state survives that
NGC6278 did not exercise. Two were found already, both invisible at small scale
until specifically hunted; a third is not unthinkable.

### 5.1 A separate A/B: the re-baselining question

Distinct from the above, and worth not conflating with it. The Fortran changes
necessarily pick a different DOP853 initial step, so libraries built after them
differ from libraries built before. At NGC6278 scale that shift was comparable
to changing the integrator tolerance by one decade.

To size it at full scale, run omega Cen with the pre-change binary and with the
post-change binary, both unchunked, and compare chi2 and the derived science
products. This measures how far existing omega Cen results move, which is the
question that decides whether published numbers need revisiting. It is a
science judgement, not an engineering one, and it does not block the
implementation.

## 6. Suggested order

1. Merger into the package, with its unit test (already written).
2. `write_orblib_dot_in` parameterisation; assert byte-identical `.in` output
   for the default arguments before going further.
3. Chunked script writer plus dispatch; NGC6278 regression via the existing
   `test_orblib_chunking.py`.
4. Config validation and documentation.
5. Full-scale A/B (section 5).
6. Only then, `orblib_chunks: 'auto'`, if wanted.

Step 2 is where a silent regression would be easiest to introduce and hardest
to notice, hence the explicit byte-comparison of the generated input files
before any behaviour changes.

## 7. Out of scope

- Adaptive chunk counts (section 2.3).
- The integrator rewrite, JAX or otherwise. If it is ever taken up, note that
  it forces the same re-baselining as section 5.1, so the two decisions are
  cheaper made together than separately.
- The integrator-tolerance systematic documented in the feasibility study
  (chi2 varying by ~1300 absolute, non-monotonically, across tolerances). It is
  pre-existing, it is unaffected by chunking, and it deserves its own
  investigation -- particularly before fine structure in a 5-10 dimensional
  chi2 landscape is trusted.
