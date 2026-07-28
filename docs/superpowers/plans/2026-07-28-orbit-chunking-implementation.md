# Implementing chunked orbit integration in DYNAMITE

Status: **implemented**, apart from the two items in section 9. Created
2026-07-28.
Feasibility study and measurements: `2026-07-28-orbit-integration-parallelisation.md`.

Section 9 records what was built and what was learned that this plan did not
anticipate. The plan itself is left as written, so that where a prediction was
wrong that is visible rather than edited away.

The Fortran side was already done and committed on this branch when this plan
was written: an orbit library can be built by N processes over disjoint orbit
ranges, and the result is bit-identical to a single-process run. What remained
was exposing that through `orblib.py` and the configuration file, documenting
it, and validating it at full scale. All but the full-scale validation is now
done; see section 9.

---

## 1. Goal and acceptance criteria

**Goal.** Let a user split each orbit family across N processes with one
configuration setting, and have the result be indistinguishable from today's
single-process library.

**The claim to defend**, at every chunk count:

1. **byte-identical** `orblib_qgrid.dat`, `orblib_losvd_hist.dat` and the box
   equivalents, after merging;
2. **identical** `chi2` and `kinchi2` from the weight solver;
3. a wall-clock reduction consistent with the chunk count.

(1) is the strong claim and the one worth designing around: if the bytes match,
everything downstream follows for free.

**Where it is checked.** NGC6278 is the development gate (section 5). It runs
in ~20 s per library on a laptop, has stored chi2 references, and already has a
harness. Omega Cen cannot serve this purpose -- one library is 9.28 CPU-h and
needs the cluster -- so it is a post-hoc confirmation (section 6), run once the
implementation is already believed correct. Nothing in development should block
on it.

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
validation protocol in section 6 (use separate output trees).

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
  makes the A/B in section 6 meaningful.

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

> **Superseded.** `auto` was implemented after all, and turned out to be cheap:
> `ModelInnerIterator.run_iteration` already computes the number of distinct
> orbit libraries an iteration will build, before any process pool starts. See
> section 9.1.

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

## 5. Development gate: NGC6278

Everything below runs on a laptop in under a minute per library, so it can be
run after every step in section 7.

### 5.1 The primary check is an invariant, not a reference value

`chunked == unchunked, same binary` is the property that matters, and it needs
no reference file at all:

- merged `_qgrid.dat` and `_losvd_hist.dat` **byte-identical** to a
  single-process run, for both families;
- `chi2`/`kinchi2` identical;
- chunk counts **1, 2, 4, 7, 8**. Seven is there deliberately: 480 does not
  divide by 7, so it exercises the remainder distribution, which is the most
  likely place for an off-by-one to hide.

Also assert the error paths: `orblib_chunks` greater than the orbit count is
clamped with a warning, and `orblib_chunks: 1` takes the old code path.

`dev_tests/test_orblib_chunking.py` already does the comparison; it needs
extending to drive the config setting rather than the executable directly.

### 5.2 `test_nnls.py` must run, and there is a trap

`test_nnls.py` does a bare `import dynamite`, so **it tests whatever is
installed, not the working tree.** Run from the repo without care, it silently
exercises site-packages and reports success against code you did not change.
It must be run as::

    PYTHONPATH=$HOME/research/dynamite python test_nnls.py user_test_config_ml.yaml

or against an editable install. Mixing the two is worse than either: installed
Python writes a 21-line `parameters_pot.in` while the repo's Fortran expects 22
(it added a Hubble constant), so the combination dies in `orbitstart` with an
EOF in `iniparam_f.f90`. Anyone bisecting this behaviour should know that up
front.

### 5.3 The tracked references are already stale

Measured with the stock 6/5/4 config, `chi2` for the three models against
`data/chi2_compare_ml_654.dat`:

| combination | model 0 | model 1 | model 2 |
|---|---|---|---|
| installed dynamite + installed binaries | -0.31% | -2.39% | -0.92% |
| **repo HEAD + pre-change binary** | **+3.62%** | **+8.34%** | **+4.37%** |
| repo HEAD + chunking binary (vs the row above) | -3.94% | -12.33% | -4.77% |

Two separate things are going on, and they must not be conflated:

1. **The repo already diverges from the reference**, before any work on this
   branch, by up to +8.3%. The installed release nearly reproduces it. The
   likely cause is the potential/dark-halo change that added `H` to
   `parameters_pot.in`. This is not ours to fix here, but it means
   `chi2_compare_ml_654.dat` cannot be used as a pass/fail gate as things
   stand.
2. **The chunking commit shifts results further**, as expected from the DOP853
   initial-step change.

Note the 6/5/4 config integrates only 120 orbits, so its Monte Carlo noise is
much larger than the 480-orbit tree used elsewhere; that is why these
percentages are bigger than the ~1-2% measured at 10/8/6. It is a
smoke-test-grade configuration, not a precision one.

**Action.** Regenerate the references from repo HEAD with
`create_comparison_data`, commit the old and new values side by side with a
note on what changed and why, and keep the invariant in 5.1 as the actual
regression gate. `test_nnls.py` is a visual test -- it plots calculated against
reference and asserts nothing -- so its role here is "the pipeline runs
end-to-end and the numbers are sane and correctly ordered", not "the numbers
match to the digit".

`dev_tests/test_orbit_losvds.py` compares against a tracked
`data/comparison_losvd.npz` and is affected the same way; regenerate it in the
same pass.

## 6. Full-scale A/B validation, after the fact

Omega Cen cannot gate development -- one library is 9.28 CPU-h and needs the
cluster -- so this runs once the NGC6278 gate in section 5 is green, to confirm
at production scale what has already been demonstrated at small scale.

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

### 6.1 A separate A/B: the re-baselining question

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

## 7. Suggested order

1. Merger into the package, with its unit test (already written).
2. `write_orblib_dot_in` parameterisation; assert byte-identical `.in` output
   for the default arguments before going further.
3. Chunked script writer plus dispatch; the section 5.1 invariant across chunk
   counts 1, 2, 4, 7, 8.
4. Config validation and documentation.
5. Regenerate `chi2_compare_ml_654.dat` and `comparison_losvd.npz` from repo
   HEAD (section 5.3), recording old and new values together.
6. `test_nnls.py` and `test_orbit_losvds.py` end to end, **with `PYTHONPATH`
   set** (section 5.2), unchunked and chunked. The two must agree with each
   other exactly; agreement with the regenerated references is then automatic.
7. Full-scale A/B (section 6), on the cluster.
8. Only then, `orblib_chunks: 'auto'`, if wanted.

Step 2 is where a silent regression would be easiest to introduce and hardest
to notice, hence the explicit byte-comparison of the generated input files
before any behaviour changes.

Step 6 is the one to be pedantic about: run it the wrong way and it passes
against the installed package while telling you nothing about the working
tree.

## 8. Out of scope

- Adaptive chunk counts (section 2.3).
- The integrator rewrite, JAX or otherwise. If it is ever taken up, note that
  it forces the same re-baselining as section 6.1, so the two decisions are
  cheaper made together than separately.
- The integrator-tolerance systematic documented in the feasibility study
  (chi2 varying by ~1300 absolute, non-monotonically, across tolerances). It is
  pre-existing, it is unaffected by chunking, and it deserves its own
  investigation -- particularly before fine structure in a 5-10 dimensional
  chi2 landscape is trusted.

---

# 9. Outcome

## 9.1 Built

- `dynamite/orblib_chunks.py` -- merges chunked libraries; byte-identical to a
  single-process run.
- `orblib_chunks` in `multiprocessing_settings`, plus `total_cores` and
  `orblib_chunks: auto`. Default 1, taking the pre-existing code path untouched.
- `LegacyOrbitLibrary.get_orbit_library_chunked` and the chunked script writer;
  `write_orblib_dot_in` gained `start`/`number`/`tag`.
- `ModelInnerIterator.resolve_orblib_chunks`, resolving the count per iteration
  from the distinct-orbit-library count, which `run_iteration` already computes
  before any pool starts.
- `dev_tests/test_orblib_chunking.py`, driving the config setting across chunk
  counts 1, 2, 4, 7, 8.
- Documented in `configuration.rst`, the API docs and the changelog.

## 9.2 What the plan got wrong

**Section 5.1's invariant was the right gate**, and it earned its place: the
first working implementation passed the chunk-count matrix and still returned
NaN for two of three models under `test_nnls`. The reason was concurrency.
Models sharing an orbit library are evaluated concurrently, so several
processes integrate into the same `datfil`. The unchunked path survives that
because every process writes the same filenames with identical content and the
last writer wins; chunked integration has more steps, and *every* shared
intermediate between them was a collision -- the chunk files, the script itself,
the merged `.dat`, and a `.bz2` that `bzip2` had created but not finished
writing.

Everything intermediate is now per-process, and only the final `.bz2` has a
shared name, appearing atomically by renaming a staging file. Three separate
attempts were needed before enumerating every filename in the path rather than
fixing them one at a time; that enumeration is what found the last two.

**This also revealed something worth its own attention:** models sharing an
orbit library each integrate it redundantly. Three models, three full
integrations. That is pre-existing, wasteful, and independent of chunking.

**The `H` attribution in section 5.3 was wrong.** The repo's deviation from
`chi2_compare_ml_654.dat` is not caused by the Hubble parameter: at the default
`H = 70` the new code writes `70 * 1e-6`, which is bit-identical to the literal
`7.0d-5` the old Fortran hardcoded, so `rho_crit` is unchanged. The same is
true of the fork-local polar-grid revert `c7eb8f5`, whose config sets 10/6/6 --
exactly the values previously hardcoded. The real cause is that the reference
was last regenerated in **#291**, before #513 (GH systematic errors now always
applied, and the test config sets `GH_sys_err`), #515 and #517 (new projected
and intrinsic mass calculations) and #442. Accumulated intended change, not a
regression, which makes regenerating it routine.

## 9.3 A memory regression found along the way

The vectorised LOSVD read shipped earlier on this branch peaked at **4.9x** the
histograms it returns, where the loop it replaced peaked at 1.01x. Its `chunk`
parameter counted (orbit, aperture) pairs while the temporaries scale with
pairs x values-per-pair, so it bounded nothing for a small library and bounded
at roughly 20x the intent for a large one. Batches are now sized by values;
peak is 2.13x on NGC6278 and the parser's own overhead is 0.41x at omega Cen
scale, with the 24x speedup unchanged.

This matters here because peak memory during the read, not the solver, is what
limits concurrent weight solves. Production monitoring shows single processes
reaching **99.6 GB** and the system reaching 636 GB.

## 9.4 Still to do

1. **Regenerate `chi2_compare_ml_654.dat` and `comparison_losvd.npz`** from repo
   HEAD, after confirming with the owners of #513/#515/#517 that the shift is
   expected. Until then `test_nnls` cannot be a pass/fail gate.
2. **Full-scale omega Cen A/B** (section 6), on the cluster.

Also worth doing, and independent of this work: re-measure the read's peak
memory on a Linux node to set `ncpus_weights` from real numbers rather than a
macOS extrapolation, and measure delta-chi2 between adjacent grid points to see
whether the integration systematic in the feasibility study matters at the
spacing an actual search uses.
