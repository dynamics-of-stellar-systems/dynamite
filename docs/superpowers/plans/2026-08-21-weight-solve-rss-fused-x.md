# Weight-Solve RSS Reduction: Fused-X Assembly, Streamed Reads, Concurrency — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this
> plan task-by-task (subagent dispatch is unreliable in this environment — run tasks inline).
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut per-worker peak RSS of a production omega-Cen adelie weight solve from ~2.5× the
NNLS matrix to ~1.2× (float64) / ~0.7× (float32) so that 4-6+ solves fit concurrently on the
1.4 TB node, with bit-identical weights and round-off-level chi2 agreement against the recorded
baseline.

**Architecture:** Three phases per the approved spec
(`docs/superpowers/specs/2026-08-21-weight-solve-rss-concurrency-design.md`):
(1) baseline RSS profiling + Linux reclaim probe on the real xeast orblib;
(2) structural changes — assemble adelie's augmented matrix X directly (A never exists) and
stream histogram reads set-by-set;
(3) stack float32, re-profile, publish the concurrency table.
Every numerical change is gated: synthetic bitwise tests first, then real-library digests, then a
full-scale end-to-end rerun compared against `PM_grid/NGC5139_adelie_xeast_output/all_models.ecsv`.

**Tech Stack:** Python 3.12, numpy>=1.26, adelie 1.1.52, scipy; dynamite repo at
`/nexus/posix0/MIA-astro-env/nneum/pesmith/dynamite`, branch `feature/orblib-performance`.
Line numbers below refer to HEAD `7b26554` (+ spec commit `1680711`).

## Global Constraints

- **PYTHONPATH trap** (documented in commit dffc9e2): installed dynamite in ENV is *not*
  editable. Every test/driver run must use
  `PYTHONPATH=/nexus/posix0/MIA-astro-env/nneum/pesmith/dynamite` so `import dynamite` resolves
  to the working tree. Verify once per session:
  `/nexus/posix0/MIA-astro-env/nneum/pesmith/ENV/bin/python -c "import dynamite; print(dynamite.__file__)"`
  must print the repo path.
- **BLAS stability**: prefix long-running commands on this KVM host with `OPENBLAS_CORETYPE=Haswell`
  (AVX-512 segfault history, see PM_grid notes).
- **Scratch space**: never write large files under `/tmp` (25 GB root disk). All scratch output
  trees live under `/nexus/posix0/MIA-astro-env/nneum/pesmith/PM_grid/`.
- **Baseline ground truth** (never modify the originals): chi2_tot = `2770835.03357815`,
  chi2_kin = `335126.55535470234`, weights at
  `PM_grid/NGC5139_adelie_xeast_output/models/orblib_000_000/ml02.60/orbit_weights.ecsv`.
- **Bit-exactness gates**: X, col_norm, y bitwise vs reference construction; returned weights
  bitwise vs baseline; chi2_tot/chi2_kin agree to ≤1e-11 relative (the final gemv is replaced by
  a residual identity — same math, different rounding; quantify the measured digits in notes).
  This refines spec gate 4 and is the agreed reading of "exact reproduction at float64".
- **scipy/cvxopt paths untouched**: `construct_nnls_matrix_and_rhs` keeps its signature and its
  role for those solvers. The only permitted edit is extracting one shared helper (Task 6).
- **Never remove the `np.shares_memory` asserts** on block writes (weight_solvers.py:875).
- **RSS-measuring runs execute sequentially** — concurrent runs would corrupt each other's
  timelines.
- Commit style follows branch history: `type(scope): summary`.

---

### Task 1: RSS tooling — reclaim probe + sampler

**Files:**
- Create: `dev_tests/_rss_probe.py`

**Interfaces:**
- Produces: `reclaim_probe(sizes_gib)` (prints a table), `class RssSampler(interval)` context
  manager with mutable `.phase` attribute writing CSV rows
  `epoch_s,phase,vmrss_kb,vmhwm_kb,anon_kb,pss_kb`. Task 2 imports both.

- [x] **Step 1: Write `_rss_probe.py`**

```python
"""Heap-reclaim probe + RSS sampler for the weight-solve memory work.

Two independent tools:

- ``reclaim_probe``: allocate arrays of given sizes, free them, report how much
  VmRSS the OS took back after del + gc.collect() (+ malloc_trim(0)). Answers
  whether Linux glibc reclaims huge numpy allocations here; the earlier macOS
  investigation measured ~0% reclaim at every size.
- ``RssSampler``: daemon thread sampling VmRSS/VmHWM (/proc/self/status) and
  Anonymous/PSS (/proc/self/smaps_rollup) into a CSV tagged with the current
  phase name.

Run standalone:
    python dev_tests/_rss_probe.py --sizes 1 8 64   # GiB
    python dev_tests/_rss_probe.py --selfcheck
"""
import argparse
import ctypes
import gc
import threading
import time

import numpy as np


def _read_status():
    out = {}
    with open('/proc/self/status') as fh:
        for line in fh:
            if line.startswith(('VmRSS:', 'VmHWM:')):
                key, val = line.split(':', 1)
                out[key] = int(val.strip().split()[0])          # kB
    return out


def _read_smaps():
    out = {}
    try:
        with open('/proc/self/smaps_rollup') as fh:
            for line in fh:
                if line.startswith(('Anonymous:', 'Pss:')):
                    key, val = line.split(':', 1)
                    out[key] = int(val.strip().split()[0])
    except OSError:
        pass
    return out


class RssSampler(threading.Thread):
    def __init__(self, path, interval=1.0):
        super().__init__(daemon=True)
        self.path = path
        self.interval = interval
        self.phase = 'init'
        self._stop_evt = threading.Event()

    def __enter__(self):
        self._fh = open(self.path, 'w', buffering=1)
        self._fh.write('epoch_s,phase,vmrss_kb,vmhwm_kb,anon_kb,pss_kb\n')
        self.start()
        return self

    def __exit__(self, *exc):
        self._stop_evt.set()
        self.join(timeout=2 * self.interval + 1)
        self._fh.close()

    def _row(self):
        st = _read_status()
        sm = _read_smaps()
        self._fh.write(
            f'{time.monotonic():.1f},{self.phase},'
            f"{st.get('VmRSS', -1)},{st.get('VmHWM', -1)},"
            f"{sm.get('Anonymous', -1)},{sm.get('Pss', -1)}\n")

    def run(self):
        while not self._stop_evt.is_set():
            self._row()
            self._stop_evt.wait(self.interval)


def _malloc_trim():
    try:
        return bool(ctypes.CDLL('libc.so.6').malloc_trim(0))
    except OSError:
        return False


def reclaim_probe(sizes_gib):
    gib = 2**30
    hdr = f'{"GiB":>6} {"rss@alloc":>10} {"after_del":>10} ' \
          f'{"after_gc":>10} {"after_trim":>10} reclaimed%'
    print(hdr)
    for s in sizes_gib:
        base = _read_status()['VmRSS']
        arr = np.empty(int(s * gib // 8), dtype=np.float64)
        arr[:: max(len(arr) // 1024, 1)] = 1.0        # touch pages sparsely
        rss_alloc = _read_status()['VmRSS']
        del arr
        rss_del = _read_status()['VmRSS']
        gc.collect()
        rss_gc = _read_status()['VmRSS']
        _malloc_trim()
        rss_trim = _read_status()['VmRSS']
        gained = rss_alloc - rss_trim                  # kB given back
        pct = 100.0 * gained / max(rss_alloc - base, 1)
        print(f'{s:>6} {rss_alloc/2**20:>9.2f}G {rss_del/2**20:>9.2f}G '
              f'{rss_gc/2**20:>9.2f}G {rss_trim/2**20:>9.2f}G {pct:>8.1f}%')


def _selfcheck():
    import os
    with RssSampler('/tmp/_rss_selfcheck.csv', interval=0.2) as samp:
        samp.phase = 'alloc'
        junk = np.empty(2**28, dtype=np.uint8); junk[:] = 1
        time.sleep(0.6)
        samp.phase = 'freed'
        del junk; gc.collect()
        time.sleep(0.6)
    lines = open('/tmp/_rss_selfcheck.csv').read().strip().splitlines()
    assert len(lines) > 4, lines
    phases = {ln.split(',')[1] for ln in lines[1:]}
    assert {'init', 'alloc', 'freed'} <= phases, phases
    print(f'selfcheck OK: {len(lines)-1} rows, phases={sorted(phases)}')
    os.remove('/tmp/_rss_selfcheck.csv')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--sizes', type=float, nargs='*', default=[1, 8, 64])
    ap.add_argument('--selfcheck', action='store_true')
    args = ap.parse_args()
    if args.selfcheck:
        _selfcheck()
    else:
        reclaim_probe(args.sizes)
```

(The selfcheck writes under `/tmp` deliberately — it is kilobytes.)

- [x] **Step 2: Run self-check**

Run: `cd /nexus/posix0/MIA-astro-env/nneum/pesmith/dynamite && PYTHONPATH=. /nexus/posix0/MIA-astro-env/nneum/pesmith/ENV/bin/python dev_tests/_rss_probe.py --selfcheck`
Expected: `selfcheck OK: N rows, phases=['alloc', 'freed', 'init']`

- [x] **Step 3: Run the real reclaim probe**

Run: `PYTHONPATH=. .../ENV/bin/python dev_tests/_rss_probe.py --sizes 1 8 64 128`
Record the table verbatim — Task 3 pastes it into the notes.

- [x] **Step 4: Commit**

```bash
git add dev_tests/_rss_probe.py
git commit -m "feat(dev_tests): RSS sampler and heap-reclaim probe"
```

### Task 2: Real-solve RSS profile driver

**Files:**
- Create: `dev_tests/_real_solve_rss_profile.py`
- Create (data, not committed): `PM_grid/NGC5139_profile_scratch_output/`,
  `PM_grid/NGC5139_config_adelie_xeast_profile.yaml`

**Interfaces:**
- Consumes: `RssSampler`, `reclaim_probe` (Task 1).
- Produces: `<csv>` timeline + `<csv-stem>_summary.json` with per-phase peak VmRSS/VmHWM;
  flags `--alm-iters N` (cap; memory observation only, never numerics), `--dtype {float64,float32}`,
  `--tag`.

- [x] **Step 1: Set up the scratch output tree** (baseline tree stays pristine)

```bash
cd /nexus/posix0/MIA-astro-env/nneum/pesmith/PM_grid
df -h /nexus/posix0     # confirm ample free space (>200 GB)
cp -a NGC5139_adelie_xeast_output NGC5139_profile_scratch_output   # real copy (~13 GB), minutes
cp NGC5139_config_adelie_xeast.yaml NGC5139_config_adelie_xeast_profile.yaml
grep -n "output_directory" NGC5139_config_adelie_xeast_profile.yaml
# then edit every io_settings/output_directory occurrence:
sed -i 's#NGC5139_adelie_xeast_output#NGC5139_profile_scratch_output#g' \
    NGC5139_config_adelie_xeast_profile.yaml
```

- [x] **Step 2: Write the driver**

```python
"""Profile the RSS timeline of a production-shape NNLS.solve(), phase by phase.

Drives the real pipeline exactly like solve() does (orblib read -> matrix ->
ALM -> final chi2), wrapping each stage so a background RssSampler tags the
timeline. Also counts ALM iterations via the bvls entry point.

The --config yaml must point at a SCRATCH copy of the baseline output tree so
solve()'s weight-file write cannot touch recorded evidence. --alm-iters caps
the multiplier loop for fast memory observation only (per-iteration memory is
flat); never use capped runs for numerical comparisons.

Run (from the dynamite repo root):
    OPENBLAS_CORETYPE=Haswell OMP_NUM_THREADS=24 \
    PYTHONPATH=$PWD \
    ENV/bin/python dev_tests/_real_solve_rss_profile.py \
      /nexus/.../PM_grid/NGC5139_config_adelie_xeast_profile.yaml \
      --tag base_f64_cap30 --alm-iters 30
"""
import argparse
import json
import resource
import sys
import time

import numpy as np

import dynamite as dyn
from dynamite import weight_solvers as ws
from _rss_probe import RssSampler      # script dir is on sys.path when run by path


def rss_gib():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == 'darwin' else r / 2**20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('config')
    ap.add_argument('--tag', required=True)
    ap.add_argument('--alm-iters', type=int, default=None)
    ap.add_argument('--dtype', default=None, choices=['float32', 'float64'])
    args = ap.parse_args()

    csv_path = (f'/nexus/posix0/MIA-astro-env/nneum/pesmith/PM_grid/'
                f'rss_{args.tag}.csv')

    c = dyn.config_reader.Configuration(args.config, reset_logging=True)
    model = c.all_models.get_model_from_row(0)
    solver = dyn.weight_solvers.NNLS(config=c, model=model)
    solver.nnls_solver = 'adelie'
    if args.dtype == 'float32':
        solver.nnls_dtype = np.float32
    if args.alm_iters is not None:
        solver.adelie_alm_iters = args.alm_iters

    orblib = model.get_orblib()
    phases = {}
    cur = {'name': 'startup'}

    def mark(name):
        phases[name] = {'t': time.perf_counter(), 'peak_rss_gib': rss_gib()}
        print(f"[phase] {name}: t={phases[name]['t']:.1f}s "
              f"peak_rss={phases[name]['peak_rss_gib']:.1f}GiB", flush=True)
        cur['name'] = name

    # wrap instance methods so the real call graph is what gets attributed
    def wrap(obj, name, label):
        orig = getattr(obj, name)
        def wrapper(*a, **k):
            mark(f'{label}:enter')
            try:
                return orig(*a, **k)
            finally:
                mark(f'{label}:exit')
        setattr(obj, name, wrapper)

    wrap(orblib, 'read_vel_histograms', 'orblib_read')
    ctor = ('construct_adelie_matrix_and_rhs'
            if hasattr(solver, 'construct_adelie_matrix_and_rhs')
            else 'construct_nnls_matrix_and_rhs')
    wrap(solver, ctor, 'matrix_build')
    wrap(solver, 'solve_adelie_alm', 'alm_solve')
    wrap(solver, 'kkt_violation' if hasattr(type(solver), 'kkt_violation_augmented')
         is False else 'kkt_violation_augmented', 'kkt')

    # count ALM iterations without touching library code
    alm_calls = {'n': 0}
    orig_bvls = ws._adelie_solver.bvls
    def counting_bvls(*a, **k):
        alm_calls['n'] += 1
        if alm_calls['n'] % 10 == 1:
            mark(f'alm_iter_{alm_calls["n"]}')
        return orig_bvls(*a, **k)
    ws._adelie_solver.bvls = counting_bvls

    t0 = time.perf_counter()
    with RssSampler(csv_path) as samp:
        samp.phase = cur['name']
        weights, chi2_tot, chi2_kin, chi2_kinmap = solver.solve(
            orblib, ignore_existing_weights=True)
        samp.phase = 'done'
    wall = time.perf_counter() - t0
    ws._adelie_solver.bvls = orig_bvls

    summary = {
        'tag': args.tag, 'wall_s': wall, 'peak_rss_gib': rss_gib(),
        'alm_iters_run': alm_calls['n'],
        'chi2_tot': None if np.isnan(chi2_tot) else float(chi2_tot),
        'chi2_kin': None if np.isnan(chi2_kin) else float(chi2_kin),
        'phases': phases,
    }
    with open(csv_path.replace('.csv', '_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    np.save(f'/nexus/posix0/MIA-astro-env/nneum/pesmith/PM_grid/'
            f'weights_{args.tag}.npy', weights)


if __name__ == '__main__':
    main()
```

(Note the kkt wrap line: after Task 5 lands the method is `kkt_violation_augmented`; before that
it falls back to `kkt_violation`. The conditional expression above selects whichever exists.)

- [x] **Step 3: Smoke-test the driver cheaply**

Run: `OPENBLAS_CORETYPE=Haswell PYTHONPATH=$PWD ENV/bin/python dev_tests/_real_solve_rss_profile.py <config> --tag smoke --alm-iters 1`
Expected: completes; summary JSON exists; `alm_iters_run` == 1. Delete `weights_smoke.npy` and the
smoke csv afterwards.

- [x] **Step 4: Commit**

```bash
git add dev_tests/_real_solve_rss_profile.py
git commit -m "feat(dev_tests): phase-attributed RSS profile driver for NNLS.solve()"
```

### Task 3: Baseline measurements → notes doc

**Files:**
- Run only: driver from Task 2, probe from Task 1
- Create: `dev_notes/weight_solve_rss_profile.md`

- [x] **Step 1: Run the capped baseline profile (sequential, nothing else running)**

```bash
cd /nexus/posix0/MIA-astro-env/nneum/pesmith/dynamite
nohup env OPENBLAS_CORETYPE=Haswell OMP_NUM_THREADS=24 PYTHONPATH=$PWD \
  /nexus/posix0/MIA-astro-env/nneum/pesmith/ENV/bin/python \
  dev_tests/_real_solve_rss_profile.py \
  /nexus/posix0/MIA-astro-env/nneum/pesmith/PM_grid/NGC5139_config_adelie_xeast_profile.yaml \
  --tag base_f64_cap30 --alm-iters 30 > /nexus/.../PM_grid/rss_base_f64_cap30.log 2>&1 &
```

Monitor: `tail -f PM_grid/rss_base_f64_cap30.log`. Expected duration ≈ 1 h.

- [x] **Step 2: Write `dev_notes/weight_solve_rss_profile.md`** containing:
  1. Measured bytes of `orblib.vel_histograms[i].y`, `intrinsic_masses`, `projected_masses` per
     set (from a tiny snippet iterating the attributes after read; include it in the notes).
  2. The Task-1 reclaim-probe table.
  3. Per-phase peaks from the summary JSON + where global peak occurs
     (assembly overlap vs X build vs ALM plateau = "resting RSS").
  4. Transient during combine/mirror at this scale (orblib_read phase delta).
  5. Derived per-worker budget: X bytes at float64/float32, projected concurrency at
     0.85×1416 GB safe budget, and the streaming go/no-go number (co-residency GiB saved).

- [x] **Step 3: Commit**

```bash
git add dev_notes/weight_solve_rss_profile.md
git commit -m "docs(notes): baseline weight-solve RSS profile on the cluster node"
```

### Task 4: Extract `_build_augmented_X` (pure refactor)

**Files:**
- Modify: `dynamite/weight_solvers.py` (solve_adelie_alm body, lines 1119-1136)
- Test: `dev_tests/test_augmented_build.py`

**Interfaces:**
- Produces: `NNLS._build_augmented_X(A_rest, b_rest, sqrt_mu, dtype)` staticmethod returning
  `(X, col_norm, y)` — bitwise what lines 1130-1136 do today.

- [x] **Step 1: Write failing-ish test first** (`test_augmented_build.py`):

```python
"""Bit-identity of the extracted augmented-matrix builder vs the inline code
it replaces in solve_adelie_alm. Run: PYTHONPATH=. python dev_tests/test_augmented_build.py
"""
import numpy as np
from dynamite.weight_solvers import NNLS


def _reference(A_rest, b_rest, sqrt_mu, dtype):
    X = np.empty((A_rest.shape[0] + 1, A_rest.shape[1]), dtype=dtype, order='F')
    X[0, :] = sqrt_mu
    X[1:, :] = A_rest
    col_norm = np.linalg.norm(X, axis=0)
    col_norm[col_norm == 0] = 1.0
    X /= col_norm
    y = np.concatenate([[0.0], b_rest]).astype(dtype)
    return X, col_norm, y


def test_bitwise():
    rng = np.random.default_rng(42)
    for dtype in (np.float64, np.float32):
        A_rest = rng.standard_normal((37, 11)).astype(dtype)
        A_rest[:, 3] = 0.0                      # exercise the zero-column guard
        b_rest = rng.standard_normal(37).astype(dtype)
        sqrt_mu = np.sqrt(1e7)
        got = NNLS._build_augmented_X(A_rest, b_rest, sqrt_mu, dtype)
        ref = _reference(A_rest, b_rest, sqrt_mu, dtype)
        for g, r in zip(got, ref):
            assert g.dtype == r.dtype and g.shape == r.shape
            assert np.array_equal(g, r), dtype


if __name__ == '__main__':
    test_bitwise()
    print('test_augmented_build OK')
```

Run it: expected FAIL with AttributeError (method missing).

- [x] **Step 2: Implement** — move lines 1122-1136 into:

```python
    @staticmethod
    def _build_augmented_X(A_rest, b_rest, sqrt_mu, dtype):
        """Build adelie's augmented, column-scaled matrix from A's body.

        Bitwise identical to the inline construction this was extracted from
        (weight_solvers.py solve_adelie_alm, pre-2026-08-21). Kept as a method
        so the fused constructor (construct_adelie_matrix_and_rhs) shares the
        exact same finishing steps."""
        n_orbs = A_rest.shape[1]
        # Filled in place into an F-ordered buffer: np.vstack + `X / col_norm`
        # + np.asfortranarray would each allocate a full copy of X, so the
        # naive version peaks at 4x the matrix (~500 GiB for omega Cen).
        X = np.empty((A_rest.shape[0] + 1, n_orbs), dtype=dtype, order='F')
        X[0, :] = sqrt_mu
        X[1:, :] = A_rest
        col_norm = np.linalg.norm(X, axis=0)
        col_norm[col_norm == 0] = 1.0
        X /= col_norm
        y = np.concatenate([[0.0], b_rest]).astype(dtype)
        return X, col_norm, y
```

and replace the moved block in solve_adelie_alm with
`X, col_norm, y = self._build_augmented_X(A_rest, b_rest, sqrt_mu, dtype)`
(`A_rest, b_rest = A[1:, :], b[1:]` stays).

- [x] **Step 3: Run** new test + `python dev_tests/test_alm_chi2_from_resid.py` → both PASS.

- [x] **Step 4: Commit** — `refactor(weight-solvers): extract _build_augmented_X from the ALM loop`

### Task 5: Surrogate chi2/kkt helpers

**Files:**
- Modify: `dynamite/weight_solvers.py` (module-level function + staticmethod beside kkt_violation)
- Test: `dev_tests/test_surrogate_chi2_kkt.py`

**Interfaces:**
- Produces:
  - `chi2_vector_from_residuals(resid_full, row0_sq)` module function → `(n_rows,)` array.
  - `NNLS.kkt_violation_augmented(row0_vec, b0, X_scaled, col_norm, resid_full, weights, mu)`
    → `(scaled, raw)`.

- [x] **Step 1: Write tests** comparing surrogate vs stock `kkt_violation(A, b, w)`:

```python
"""Surrogate chi2-vector and KKT identities vs direct computation.

The surrogates are algebraically identical to A-based forms but re-order the
floating-point ops, so they are tested to rtol 1e-10, NOT bitwise. The scaled
KKT value must stay inside [0, 1].
Run: PYTHONPATH=. python dev_tests/test_surrogate_chi2_kkt.py
"""
import numpy as np
from scipy.optimize import nnls
from dynamite.weight_solvers import NNLS, chi2_vector_from_residuals


def _problem(rng, n_rows, n_orbs, dtype=np.float64):
    A = rng.standard_normal((n_rows, n_orbs)).astype(dtype)
    A[0, :] = 1e8                       # mimic the total-mass row scale
    w_true = rng.random(n_orbs) ** 3    # mostly-near-zero weights
    w_true[: 3] = 0.0                   # force some exact zeros
    b = A @ w_true.astype(dtype)
    return A, b


def test_kkt_matches_stock():
    rng = np.random.default_rng(7)
    A, b = _problem(rng, 300, 40)
    mu = 1e7
    X, col_norm, _y = NNLS._build_augmented_X(A[1:], b[1:], np.sqrt(mu),
                                              np.float64)
    row0_vec = np.full(A.shape[1], 1.0) / max(abs(b[0]), 1.0)  # arbitrary row0
    w, _ = nnls(A, b)
    resid_full = A @ w - b              # plain residual, aligned to A rows
    got = NNLS.kkt_violation_augmented(row0_vec, b[0], X, col_norm,
                                       resid_full, w, mu)
    ref = NNLS.kkt_violation(A, b, w)
    assert 0.0 <= got[0] <= 1.0
    assert np.isclose(got[0], ref[0], rtol=1e-10), (got, ref)
    assert np.isclose(got[1], ref[1], rtol=1e-8), (got, ref)


def test_exact_fit_returns_zero_scaled():
    rng = np.random.default_rng(11)
    n_orbs, n_mass = 12, 30
    A = rng.standard_normal((n_mass, n_orbs))
    w = np.zeros(n_orbs); w[:5] = rng.random(5) + 0.1
    b = A @ w                            # exactly attainable
    mu = 1e7
    X, col_norm, _ = NNLS._build_augmented_X(A[1:], b[1:], np.sqrt(mu),
                                             np.float64)
    row0_vec = A[0]
    resid_full = A @ w - b               # ~0 everywhere
    scaled, raw = NNLS.kkt_violation_augmented(row0_vec, b[0], X, col_norm,
                                               resid_full, w, mu)
    assert raw <= 1e-8 and scaled <= 1e-8


def test_chi2_vector_identity():
    rng = np.random.default_rng(13)
    A, b = _problem(rng, 120, 25)
    w, _ = nnls(A, b)
    resid_full = A @ w - b
    got = chi2_vector_from_residuals(resid_full, float(resid_full[0]))
    ref = (A @ w - b) ** 2
    assert np.allclose(got, ref, rtol=0, atol=0)   # pure reshape/square: bitwise
    assert got.shape == ref.shape


if __name__ == '__main__':
    test_kkt_matches_stock()
    test_exact_fit_returns_zero_scaled()
    test_chi2_vector_identity()
    print('test_surrogate_chi2_kkt OK')
```

Note `chi2_vector_from_residuals(resid_full, row0_sq)` returns
`np.concatenate(([row0_sq], resid_full[1:] ** 2))` — index 0 slot is supplied by the caller's
row-0 term; the test passes `resid_full[0]**2` as that term, making the comparison exact.

- [x] **Step 2: Run → FAIL (functions missing). Step 3: implement:**

```python
def chi2_vector_from_residuals(resid_full, row0_sq):
    """Per-row squared residuals of ``A @ w - b`` without materializing A.

    ``resid_full`` holds rows 1..n of that residual (from adelie's state.resid
    upcast to float64); the total-mass row's contribution arrives separately
    as ``row0_sq`` because X replaces that row. Index 0 of the result is
    ``row0_sq`` itself, keeping the indexing of the A-based chi2_vector
    (chi2_kin slices ``[n_mass_constraints:]``)."""
    resid_full = np.asarray(resid_full, dtype=np.float64).ravel()
    return np.concatenate(([row0_sq], resid_full[1:] ** 2))
```

and on NNLS (beside kkt_violation):

```python
    @staticmethod
    def kkt_violation_augmented(row0_vec, b0, X_scaled, col_norm, resid_full,
                                weights, mu):
        r"""kkt_violation expressed through the augmented matrix.

        Identities (col_norm includes the penalty row):
            A[1:, j]      = col_norm[j] * X_scaled[1:, j]
            grad[j]       = row0_vec[j]*r0 + col_norm[j]*(X_scaled[1:]^T r)[j]
            ||A_.j||^2    = row0_vec[j]^2 + col_norm[j]^2 - mu
        Two passes over X instead of three over A; no A anywhere. Same
        semantics as kkt_violation: returns (scaled in [0,1], raw).
        """
        r0 = float(row0_vec @ weights) - float(b0)
        r_rest = np.asarray(resid_full, dtype=np.float64).ravel()[1:]
        resid = np.concatenate(([r0], r_rest))
        grad = (row0_vec.astype(np.float64) * r0
                + col_norm.astype(np.float64)
                * (X_scaled[1:, :].T @ r_rest))
        viol = np.where(weights > 0, np.abs(grad), np.maximum(-grad, 0.0))
        raw = float(np.max(viol))
        col_sq = np.maximum(col_norm.astype(np.float64) ** 2 - mu
                            + row0_vec.astype(np.float64) ** 2, 0.0)
        scale = np.sqrt(col_sq) * np.linalg.norm(resid)
        if not np.any(scale > 0):
            return 0.0, raw
        scale = np.where(scale > 0, scale, np.inf)
        return float(np.max(viol / scale)), raw
```

IMPORTANT subtlety encoded in the tests: the caller supplies the *plain* residual aligned to A
rows. Inside solve_adelie_alm, `state.resid` is the AUGMENTED residual whose slot 0 is the penalty
row; the plain rows-1.. residual equals `state.resid[1:]` only if y[1:]==b_rest (true: y is built
as `[0, b_rest...]`). Document this at the call site wired in Task 7.

- [x] **Step 4: Run → PASS. Step 5: Commit** — `feat(weight-solvers): surrogate chi2/KKT identities for the fused path`

### Task 6: Fused constructor `construct_adelie_matrix_and_rhs`

**Files:**
- Modify: `dynamite/weight_solvers.py`
- Test: `dev_tests/test_fused_assembly.py`

**Interfaces:**
- Produces: `class AdelieProblem(NamedTuple)` fields `X, col_norm, y, row0_vec, b0`;
  `NNLS.construct_adelie_matrix_and_rhs(orblib) -> AdelieProblem`; shared helper
  `NNLS._prepare_kinematic_block(kins, orb_veldist, tmp, prj_mass_i)` used by BOTH constructors.

- [x] **Step 1: Extract `_prepare_kinematic_block`** from construct_nnls_matrix_and_rhs's loop
  body (lines ~829-849: everything between picking prj_mass_i and appending constraints):

```python
    def _prepare_kinematic_block(self, kins, orb_veldist, tmp, prj_mass_i):
        """Shared by both NNLS constructors: scale observed kinematics by the
        set's projected masses, apply the CR-cut edge-bin zeroing to the orbit
        histograms IN PLACE (mirrors triaxnnls_CRcut.f90), transform the orbit
        library into the observed parameterisation, apply the CRcut if enabled.
        Returns flat (obs_kins, obs_err, orb_kins)."""
        hist_dim = len(orb_veldist.y[0,...,0].shape)  # 1D or 2D vel hists
        obs_kins, obs_kins_err = tmp
        obs_kins = (obs_kins.T * prj_mass_i).T
        obs_kins_err = (obs_kins_err.T * prj_mass_i).T
        if hist_dim == 1:
            orb_veldist.y[:,0,:] = 0.
            orb_veldist.y[:,-1,:] = 0.
        orb_kins = kins.transform_orblib_to_observables(orb_veldist,
                                                        self.settings)
        if self.CRcut:
            orb_kins = self.apply_CR_cut(kins, orb_veldist, orb_kins)
        return np.ravel(obs_kins), np.ravel(obs_kins_err), orb_kins
```

Replace the corresponding block inside construct_nnls_matrix_and_rhs with calls to it (same
statement order ⇒ same floating point). Keep `idx_ap_start` bookkeeping at the call site.

- [x] **Step 2: Add AdelieProblem + fused constructor** (new method directly after
  construct_nnls_matrix_and_rhs). Core structure — mirror stock order exactly:

```python
class AdelieProblem(typing.NamedTuple):
    """Everything the adelie path needs, built without materializing A.

    X has exactly n_rows rows: the sqrt_mu penalty row REPLACES A's total-mass
    row, which survives only as the small vector row0_vec (= ones/econ[0],
    bitwise equal to A[0])."""
    X: np.ndarray          # (n_rows, n_orbs) F-order, column-scaled
    col_norm: np.ndarray   # (n_orbs,) unit-L2 norms incl. penalty row
    y: np.ndarray          # (n_rows,) ALM target; slot 0 updated per iteration
    row0_vec: np.ndarray   # (n_orbs,) == A[0] bitwise
    b0: float              # rhs[0] = total_mass/total_mass_error
```

```python
    def construct_adelie_matrix_and_rhs(self, orblib):
        """Assemble adelie's augmented matrix directly — A never exists.

        Writes the sqrt(mu) penalty row into X row 0, streams every constraint
        block straight into rows 1.., divides by econ, then finishes with the
        SAME col_norm/divide/y steps as _build_augmented_X, so X/col_norm/y are
        bit-identical to the two-step construction. Saves one full matrix
        (~125 GiB for omega Cen)."""
        dtype = self.nnls_dtype
        stars = self.system.get_unique_triaxial_visible_component()
        obs_values = [kins.get_observed_values_and_uncertainties(self.settings)
                      for kins in stars.kinematic_data]
        n_rows = self.n_mass_constraints + sum(np.size(v) for v, _ in obs_values)
        n_orbs = orblib.n_orbs
        sqrt_mu = np.sqrt(self.adelie_mu)
        con = np.zeros(n_rows, dtype=dtype)
        econ = np.zeros(n_rows, dtype=dtype)
        X = np.empty((n_rows, n_orbs), dtype=dtype, order='F')
        X[0, :] = sqrt_mu
        con[0] = self.total_mass
        econ[0] = self.total_mass_error
        idx = slice(1, 1 + self.n_intrinsic)
        con[idx] = np.ravel(self.intrinsic_masses)
        error = np.abs(np.ravel(self.intrinsic_masses
                                * self.intrinsic_mass_error))
        error[np.where(error <= 0.0)] = 1.0e-16
        econ[idx] = error
        X[idx, :] = np.reshape(orblib.intrinsic_masses, (n_orbs, -1)).T
        idx = slice(1 + self.n_intrinsic,
                    1 + self.n_intrinsic + self.n_apertures)
        con[idx] = self.projected_masses
        econ[idx] = np.abs(self.projected_masses * self.projected_mass_error)
        X[idx, :] = np.hstack(orblib.projected_masses).T
        idx_ap_start = 0
        idx_row = self.n_mass_constraints
        for (kins, orb_veldist, tmp) in zip(stars.kinematic_data,
                                            orblib.vel_histograms, obs_values):
            n_ap = kins.n_spatial_bins
            prj_mass_i = self.projected_masses[idx_ap_start:idx_ap_start + n_ap]
            idx_ap_start += n_ap
            obs_kins, obs_kins_err, orb_kins = self._prepare_kinematic_block(
                kins, orb_veldist, tmp, prj_mass_i)
            n_orb_constraints = orb_kins.size // n_orbs
            idx_row_end = idx_row + obs_kins.size
            assert n_orb_constraints == obs_kins.size, \
                f'{type(kins).__name__}: orbit library gives ' \
                f'{n_orb_constraints} constraints per orbit but the ' \
                f'kinematic data gives {obs_kins.size}'
            con[idx_row:idx_row_end] = obs_kins
            econ[idx_row:idx_row_end] = obs_kins_err
            dest = X[idx_row:idx_row_end, :].T.reshape(orb_kins.shape)
            assert np.shares_memory(dest, X), \
                'X block write got a copy, not a view - block would be ' \
                'silently left at zero'
            dest[...] = orb_kins
            idx_row = idx_row_end
        if np.any(con[econ == 0] != 0):
            raise ValueError('Weight solving fail: zero errors for nonzero '
                             'constraints!')
        rhs = np.zeros_like(con)
        np.divide(con, econ, out=rhs, where=econ != 0)
        # zero-error fixup on matrix rows 1.. BEFORE division, mirroring the
        # stock constructor (its row-0 equivalent becomes row0_vec below)
        econ_body = econ[1:]
        nz = (X[1:, :] != 0) & (econ_body == 0)[:, None]
        if np.any(nz):
            rr, oo = np.nonzero(nz)
            txt = (f'Weight solving problem in {self.direc_with_ml}: '
                   'zero errors for nonzero matrix coefficients at '
                   f'[constraint no, orbit no] = {(rr + 1, oo)}! Matrix '
                   f'value(s) there ({X[1:, :][nz]}) will be considered zero.')
            self.logger.warning(txt)
            X[1:, :][nz] = 0
        # divide rows by their errors: same elementwise op as the stock
        # transposed-view divide, just restricted to rows 1.. (row 0 of A
        # becomes row0_vec below)
        body = X[1:, :].T                          # view (n_orbs, n_rows-1)
        np.divide(body, econ_body, out=body, where=econ_body != 0)
        # finish exactly like _build_augmented_X
        col_norm = np.linalg.norm(X, axis=0)
        col_norm[col_norm == 0] = 1.0
        X /= col_norm
        y = np.concatenate([[0.0], rhs[1:]]).astype(dtype)
        row0_vec = np.full(n_orbs, 1.0, dtype=dtype) / econ[0]
        return AdelieProblem(X=X, col_norm=col_norm, y=y, row0_vec=row0_vec,
                             b0=float(rhs[0]))
```

- [x] **Step 3: Write `test_fused_assembly.py`** — a fake-system harness driving BOTH real
  constructors end-to-end and asserting bitwise equality:

Harness sketch (full file in this spirit; fake kins return `np.moveaxis(raw, -1, 1)` views so the
reshape-view assert is genuinely exercised):

```python
"""Bitwise equivalence: construct_adelie_matrix_and_rhs(orblib) ==
_build_augmented_X(*construct_nnls_matrix_and_rhs(orblib) sliced).

Uses a minimal fake system/kins/orblib (real Histogram objects, moveaxis'd
transform outputs) so the REAL methods run, including the view-write asserts.
Run: PYTHONPATH=. python dev_tests/test_fused_assembly.py
"""
import types
import numpy as np
from dynamite import kinematics as dyn_kin
from dynamite.weight_solvers import NNLS


class FakeKins:
    kind = 'gh'
    def __init__(self, shape, n_bins):
        self.shape = shape
        self.n_spatial_bins = n_bins
        self.data = ...
    def get_observed_values_and_uncertainties(self, settings):
        n_c = int(np.prod(shape))
        return ((rng stuff shaped (n_bins, ...)), (err same shape))
    def transform_orblib_to_observables(self, hist, settings):
        # non-contiguous view, like the real ones
        return np.moveaxis(hist.y[:, ..., :].transpose(...), ...)
```

(Concrete shapes: one 1D set with bins=6, velocity bins=5; one 2D PM-like set 4×5 velocity grid,
3 spatial bins. FakeOrblib carries `n_orbs`, `vel_histograms` list of real `dyn_kin.Histogram` /
`dyn_kin.Histogram2D`, `intrinsic_masses` (n_orbs, 4), `projected_masses` list. FakeNNLS built
with `types.SimpleNamespace` carrying every attribute the constructors touch — nnls_dtype,
adelie_mu, CRcut=False, settings={}, n_mass_constraints, n_intrinsic=4, n_apertures=sum(n_bins),
intrinsic/projected masses+errors (MGE side), total_mass/error, direc_with_ml='fake',
logger=module logger, system stub exposing get_unique_triaxial_visible_component() → namespace
with kinematic_data=[FakeKins...]; then bind methods:
`fused = NNLS.construct_adelie_matrix_and_rhs(fake, orblib)` etc.)

Assertions:

```python
A_ref, b_ref = NNLS.construct_nnls_matrix_and_rhs(fake, orblib)
X_ref, cn_ref, y_ref = NNLS._build_augmented_X(
    A_ref[1:], b_ref[1:], np.sqrt(mu), dtype)
prob = NNLS.construct_adelie_matrix_and_rhs(fake, orblib)
assert np.array_equal(prob.X, X_ref)
assert np.array_equal(prob.col_norm, cn_ref)
assert np.array_equal(prob.y, y_ref)
assert np.array_equal(prob.row0_vec, A_ref[0])       # bitwise A row 0
assert prob.b0 == float(b_ref[0])
```

Repeat for dtype float32. Also mutate `hist.y` edge bins between calls? No — the CR-cut zeroing
is destructive; call order above reuses the same orblib for both paths, and the FIRST
(reference) call zeroes edges already, so the second sees pre-zeroed hists — values still agree
because zeroing is idempotent. Add a comment noting this.

- [x] **Step 4: Run** `PYTHONPATH=. python dev_tests/test_fused_assembly.py` → PASS; rerun
  `dev_tests/test_nnls_matrix_assembly.py` and `dev_tests/test_augmented_build.py` → PASS
  (refactor did not drift the stock constructor).

- [x] **Step 5: Commit** — `feat(weight-solvers): fused augmented-matrix constructor (A never materialized)`

### Task 7: Rewire solve()/solve_adelie_alm onto the fused path

**Files:**
- Modify: `dynamite/weight_solvers.py` (solve_adelie_alm signature/body; solve() adelie branch +
  final chi2 block)
- Modify: `dev_tests/_real_alm_chi2_check.py` (new signature)
- Create: `dev_tests/_real_fused_check.py` (two-mode digest comparator on the real orblib)

**Interfaces:**
- `solve_adelie_alm(self, problem: AdelieProblem) -> tuple[np.ndarray, np.ndarray]`
  returning `(best_w, resid_full)` where `resid_full` is the PLAIN residual
  (`Aw-b`) aligned to A rows: slot 0 = `row0_vec@best_w - b0`, slots 1.. =
  `problem.y[1:] - problem.X[1:] @ (col_norm*best_w)`.
- `solve()` adelie branch: `problem = self.construct_adelie_matrix_and_rhs(orblib)`;
  nan fallback length `problem.X.shape[1]`; final chi2 via
  `chi2_vector_from_residuals(resid_full, row0_resid**2)` with
  `row0_resid = float(problem.row0_vec @ weights) - problem.b0`.
  scipy/cvxopt branches untouched.

- [x] **Step 1: Rewrite solve_adelie_alm** around `problem`:
  - keep `sqrt_mu`, thread caps, lower/upper, weights_arr, lam loop, best tracking unchanged;
  - `w = (np.asarray(state.beta).ravel() / col_norm).astype(np.float64)` unchanged;
  - per-iteration chi2 line becomes
    `row0 = float(problem.row0_vec @ w) - problem.b0`
    (bitwise-equal to old `float(A[0] @ w) - float(b[0])`);
  - after the loop compute `resid_full` for best_w as above (one gemv over X);
  - KKT: `self.kkt_violation_augmented(problem.row0_vec, problem.b0, X, col_norm,
    resid_full, best_w, mu)`; warning threshold block unchanged;
  - `return best_w, resid_full`;
  - update docstring: takes AdelieProblem; explain the residual identities and why the final
    chi2 moves to them (round-off-only difference vs the old gemv, quantified in notes).

- [x] **Step 2: Rewire solve()** exactly as Interfaces specifies. Keep the existing
  `if self.nnls_solver != "adelie":` normalization guard; adelie branch no longer needs A/b at
  all. Exception handler: log as today; `weights = np.full(problem.X.shape[1], np.nan)`.

- [x] **Step 3: Update `_real_alm_chi2_check.py`:** after loading A_sub/b_sub build the problem
  directly:

```python
sqrt_mu = np.sqrt(solver.adelie_mu)
dtype = A.dtype
X, col_norm, y = dyn.weight_solvers.NNLS._build_augmented_X(
    A[1:], b[1:], sqrt_mu, dtype)
problem = dyn.weight_solvers.AdelieProblem(
    X=X, col_norm=col_norm, y=y,
    row0_vec=A[0], b0=float(b[0]))
weights, resid_full = solver.solve_adelie_alm(problem)
direct = float(np.sum((A @ weights - b) ** 2))
reported = float(np.sum(dyn.weight_solvers.chi2_vector_from_residuals(
    resid_full, (float(problem.row0_vec @ weights) - problem.b0) ** 2)))
print(reported, direct, abs(reported - direct) / direct)
```

Expected: relative difference < 1e-11 at float64.

- [x] **Step 4: Create `_real_fused_check.py`** — digest comparator proving bitwise equality on
  the real orblib WITHOUT holding both matrices:

```python
"""Real-library proof that the fused constructor is bitwise the two-step one.

Two invocations (each peaks at ~1x matrix + orblib histograms, never 2x):
    ... _real_fused_check.py CONFIG --mode reference --out ref.json
    ... _real_fused_check.py CONFIG --mode fused     --out fused.json
Then any diff of the two JSONs must be empty. Digests are sha256 of contiguous
slabs (rows), so equality proves bit-identity of X; col_norm/y/row0_vec/b0 are
compared exactly via the stored hex of their bytes.
"""
```

Implementation points: drive like Task 2 (config/model/orblib/read); mode `reference` computes
`(A, b) = construct_nnls_matrix_and_rhs`, then `X, col_norm, y = _build_augmented_X(...)`,
digests slab-wise (slab = 2048 rows), records `{'shape':…, 'dtype':str(X.dtype),
'col_norm_sha':sha256(col_norm.tobytes()), 'y_sha':…, 'row0_sha':sha256(A[0].tobytes()),
'b0':b[0], 'x_slabs':[...]}`; mode `fused` builds `AdelieProblem` and digests identically.
Free A before digesting in reference mode if memory is tight (digest needs only X).

- [x] **Step 5: Run the checks**

```bash
OPENBLAS_CORETYPE=Haswell OMP_NUM_THREADS=24 PYTHONPATH=$PWD \
  ENV/bin/python dev_tests/_real_fused_check.py \
  /nexus/.../PM_grid/NGC5139_config_adelie_xeast_profile.yaml \
  --mode reference --out /nexus/.../PM_grid/fused_ref.json
# then same with --mode fused --out .../fused_new.json
diff <(python -c "import json;d=json.load(open('...'));d.pop('timing');print(json.dumps(d,indent=1,sort_keys=True))" < ref) <(same for new)
```
Expected: identical JSON (excluding timing keys). Runtime ≈ 2 × ~15 min.

Also rerun ALL unit tests: `test_augmented_build.py`, `test_surrogate_chi2_kkt.py`,
`test_fused_assembly.py`, `test_alm_chi2_from_resid.py`, `test_nnls_matrix_assembly.py` → green.

- [x] **Step 6: Commit** — `perf(weight-solvers): route the adelie path through the fused constructor`

### Task 8: Streamed per-set histogram reads

**Files:**
- Modify: `dynamite/orblib.py` (`read_vel_histograms`, `read_orbit_base`)
- Modify: `dynamite/weight_solvers.py` (streaming branch in fused constructor), config plumbing
  wherever `nnls_dtype` is read (mirror it — grounding step below)
- Modify: `docs/getting_started/configuration.rst` (one entry)
- Test: `dev_tests/test_streamed_reads.py`

**Design:** `read_vel_histograms(pops=False, kin_sets=None, skip_density=False)`.
`kin_sets=None` (default) is byte-for-byte today. With `kin_sets=[i]` only set i's histograms are
parsed (per-set aperture subsets threaded into the bulk parsers — they already accept subsets),
combined/mirrored, and stored; other `vel_histograms` entries become None.
`skip_density=True` skips the qgrid/density parse entirely and leaves
`intrinsic_masses`/`n_orbs` from the previous call (caller asserts they exist). The fused
constructor gains: when `stream_orblib_reads: true` (weight_solver_settings, default false until
validated), it (a) calls once normally for densities/masses/n_orbs, (b) loops sets calling
`read_vel_histograms(kin_sets=[i], skip_density=(i>first_1d_call))`, writing each block into X,
then `orblib.vel_histograms[i] = None`.

- [x] **Step 1: Grounding read** — orblib.py lines 1689-2095 (`read_orbit_base` + bulk parsers'
  aperture-subset arguments) and the `nnls_dtype` plumbing in weight_solvers/config_reader.
  Confirm how parser subset indices map kinematic-set → file-local apertures.

- [x] **Step 2: Write `test_streamed_reads.py`** reusing the FortranFile-writer fixtures from
  `test_hist_bulk_read.py` (import its writer helpers): build a multi-set library (≥2 1D sets +
  ≥2 2D sets, mixed order, some empty sentinels); parse fully; parse per-set; assert each
  streamed `Histogram.y` is `np.array_equal` to the corresponding full-parse slice; assert
  `skip_density` second call leaves `intrinsic_masses`/`n_orbs` untouched and skips qgrid
  (monkeypatch counter on the qgrid reader).

- [x] **Step 3: Implement** the signatures above. Careful points: aperture-index bookkeeping per
  file (sets interleave within a file — the mixed-layout lesson from test_hist_bulk_read);
  `projected_masses` computed only for parsed sets; `scale_x_values(velocity_scaling_factor)`
  applied per parsed set exactly as the full path does; mirror handling unchanged
  (`combine_and_mirror_orblibs(t, b, mirror=mirror)`).

- [x] **Step 4: Wire the streaming branch** into `construct_adelie_matrix_and_rhs` behind
  `getattr(self, 'stream_reads', False)`; document the setting in configuration.rst next to
  `nnls_dtype` ("memory setting; results identical").

- [x] **Step 5: Verify**: new test green; `test_hist_bulk_read.py`, `test_pm_hist_bulk_read.py`,
  `test_fused_assembly.py` green; rerun `_real_fused_check.py --mode fused` with
  `stream_orblib_reads: true` added to the scratch config → JSON identical to the non-streamed
  fused run (this is the real-library proof that streaming changed nothing).

- [x] **Step 6: Commit** — `feat(orblib): stream per-set histogram reads for the fused NNLS path`

### Task 9: Full-scale validation, re-profiles, concurrency table, docs

**Files:**
- Runs: driver (Task 2), checker (Task 7), probe (Task 1)
- Modify: `dev_notes/weight_solve_performance.md`, `dev_notes/weight_solve_rss_profile.md`,
  `docs/more_info/changelog.rst`

- [x] **Step 1: End-to-end float64 validation** (overnight, sequential):

```bash
cd /nexus/posix0/MIA-astro-env/nneum/pesmith/dynamite
nohup env OPENBLAS_CORETYPE=Haswell OMP_NUM_THREADS=24 PYTHONPATH=$PWD \
  ENV/bin/python dev_tests/_real_solve_rss_profile.py \
  /nexus/.../PM_grid/NGC5139_config_adelie_xeast_profile.yaml \
  --tag fused_f64_full > /nexus/.../PM_grid/rss_fused_f64_full.log 2>&1 &
```
(no `--alm-iters` ⇒ full run, ~4-5 h). Acceptance: `weights_fused_f64_full.npy` bitwise-equal to
the baseline `orbit_weights.ecsv` weights column (load ecsv with astropy, compare
`np.array_equal`); chi2_tot/kin within 1e-11 relative of
2770835.03357815 / 335126.55535470234 — record the actual digits achieved.

- [x] **Step 2: Re-profile matrix** (each capped ~1 h, strictly sequential):
  `base_f64_cap30` (already have), `fused_f64_cap30`, `fused_stream_f64_cap30`,
  `fused_f32_cap30`, `fused_stream_f32_cap30`.

- [x] **Step 3: Concurrency table** appended to `weight_solve_rss_profile.md`:
  per-mode peak GiB; recommended `ncpus_weights = floor(1200 GB / peak)` with the 0.85 safety
  factor; note OMP_NUM_THREADS sizing per worker (≈192/ncpus_weights).

- [x] **Step 4: Docs**: update `weight_solve_performance.md` ("where memory goes now" section +
  close the release-A open item as superseded-by-fusion), changelog entries (fused constructor,
  streaming reads, RSS profile tooling), and a short paragraph in the RSS notes stating the
  chi2 round-off refinement of spec gate 4 with the measured digits.

- [x] **Step 5: Commits** — `perf(weight-solvers): validate fused path at production scale`,
  `docs: concurrency table and consolidated weight-solve memory notes`.

---

## Self-review record

- Spec coverage: Phase 1 → Tasks 1-3; Phase 2a → Tasks 4-7; Phase 2b → Task 8; Phase 3 → Tasks
  8-9 (float32 flag existed already; composition validated in Task 9 step 2). Concurrency table →
  Task 9. All five spec gates mapped: gates 1→Tasks 4/6, 2→Task 5, 3→Task 7 (`_real_fused_check`),
  4→Task 9 step 1, 5→Task 9 steps 2-3.
- Placeholders: none — every task carries concrete code or exact anchors; the two "grounding
  read" steps (Task 8 step 1, Task 6 harness details) pin interfaces with concrete code while
  pointing at exact line ranges.
- Type consistency: `AdelieProblem` fields used identically in Tasks 6/7; `_build_augmented_X`
  signature consistent across Tasks 4/5/6/7; `kin_sets`/`skip_density` naming consistent in
  Task 8.
