# Plan: eliminate the read/mirror/combine copy chain in `read_vel_histograms`

Status: **implemented** (`orblib.py`: `duplicate_flip_and_interlace_orblib`
+ `combine_orblibs` replaced by `combine_and_mirror_orblibs`). Verified
against the real patched file: exact pre-change chi2 reproduced via three
independent scripts (a dedicated verification script, the `nnls_dtype`
float32/float64 validation, and the repo's pre-existing
`dev_tests/test_nnls_adelie.py`), on both `nnls_dtype` settings. The
`mirror=False` (bar/disk) and 2D-histogram (proper motions) branches are
implemented (the function branches on `orblib.y.ndim` and the `mirror` flag
explicitly) but still lack a real-orblib test the way the mirrored/1D path
now has - no locally-built orbit library exercises either case yet.

Written after measuring that a
`del`-based cleanup of `orblib.vel_histograms`/`intrinsic_masses`/
`projected_masses` does not reduce process RSS (tested at 0.36/1.2/3.5/7.7 GB,
all ~0% reclaimed - see `dynamite_production_memory.md` and the OOM
investigation this plan follows on from). The `del` approach fails because
the data is freed but not returned to the allocator/OS; this plan instead
avoids *allocating* the excess in the first place, which is unaffected by
that reclaim problem.

## Where the excess comes from

`OrbitLibrary.read_vel_histograms()` (`orblib.py:1928`) builds each
kinematic dataset's combined orbit library through three generations of full
array copies:

1. `read_orbit_base("orblib")` reads the tube orbits: one dense array per
   kinematic dataset, shape `(n_tube, n_vel, n_ap)`. This step is already
   tight - not part of the problem.
2. `duplicate_flip_and_interlace_orblib` (`orblib.py:1818`) allocates a new
   `(2*n_tube, n_vel, n_ap)` array and fills it by interlacing the tube data
   with its own velocity-reversed copy. While this runs, the original tube
   array and the new 2x array are both live.
3. `read_orbit_base("orblibbox")` reads the box orbits: another dense
   `(n_box, n_vel, n_ap)` array, while the mirrored tube array from step 2 is
   still referenced.
4. `combine_orblibs` (`orblib.py:1875`) allocates a *third* array,
   `(2*n_tube + n_box, n_vel, n_ap)`, and copies the mirrored tube data and
   the box data into it. During this copy, the mirrored tube array, the box
   array, and the new combined array are all live simultaneously.

Peak concurrent memory during this sequence is roughly 3-4x the final
combined array's size, per kinematic dataset. This is what was measured as
"a transient 2-3x peak during the read phase" earlier in the investigation.

## Proposed fix: one pre-allocated buffer, filled in place

Replace steps 2+4 with a single function that allocates the *final* combined
array once and writes both sources directly into their slices - no
intermediate "mirrored" array ever exists.

```python
def combine_and_mirror_orblibs(self, tube, box, mirror=True):
    """Build the final per-kinematic-dataset orbit library array directly,
    without materializing a separate mirrored-tube generation.

    Parameters
    ----------
    tube : Histogram or Histogram2D
        as returned by read_orbit_base("orblib")
    box : Histogram or Histogram2D
        as returned by read_orbit_base("orblibbox")
    mirror : bool
        whether to mirror+interlace the tube orbits (False for
        bar/disk systems, matching the existing is_bar_disk_system() check)

    Returns
    -------
    Histogram or Histogram2D
        same shape/content as
        combine_orblibs(duplicate_flip_and_interlace_orblib(tube), box)
        produces today, built with one allocation instead of three.
    """
    n_tube = tube.y.shape[0]
    n_box = box.y.shape[0]
    n_tube_out = 2 * n_tube if mirror else n_tube
    final_shape = (n_tube_out + n_box,) + tube.y.shape[1:]
    final = np.zeros(final_shape, dtype=tube.y.dtype)

    if mirror:
        if tube.y.ndim == 3:  # 1D histograms (losvd)
            final[0:2 * n_tube:2] = tube.y
            final[1:2 * n_tube:2] = tube.y[:, ::-1, :]
        else:  # 2D histograms (proper motions)
            final[0:2 * n_tube:2] = tube.y
            final[1:2 * n_tube:2] = tube.y[:, ::-1, ::-1, :]
    else:
        final[0:n_tube] = tube.y

    final[n_tube_out:] = box.y

    cls = dyn_kin.Histogram if tube.y.ndim == 3 else dyn_kin.Histogram2D
    return cls(xedg=tube.xedg, y=final)
```

`read_vel_histograms` changes from:

```python
tube_orblib, tube_density_3D = self.read_orbit_base("orblib", pops=pops)
if not self.system.is_bar_disk_system():
    tmp = [self.duplicate_flip_and_interlace_orblib(t) for t in tube_orblib]
    tube_orblib = tmp
    ...
box_orblib, box_density_3D = self.read_orbit_base("orblibbox", pops=pops)
self.vel_histograms = [self.combine_orblibs(t, b)
                        for t, b in zip(tube_orblib, box_orblib)]
```

to:

```python
tube_orblib, tube_density_3D = self.read_orbit_base("orblib", pops=pops)
box_orblib, box_density_3D = self.read_orbit_base("orblibbox", pops=pops)
mirror = not self.system.is_bar_disk_system()
self.vel_histograms = [
    self.combine_and_mirror_orblibs(t, b, mirror=mirror)
    for t, b in zip(tube_orblib, box_orblib)
]
if mirror:
    tube_density_3D = np.repeat(tube_density_3D, 2, axis=0)
```

Note the density_3D repeat (a separate, much smaller array - not the
memory driver) still happens the same way; only the `vel_histograms`
per-dataset arrays go through the new path.

## Expected peak reduction

Per kinematic dataset, concurrent live data during construction drops from
roughly `3-4x final_size` to `1x final_size + max(tube_size, box_size)`
(the raw tube/box reads, each needed only transiently until their slice is
written, then dead). For a dataset where tube and box are comparable in
size to the final combined array, this is close to a 2-3x reduction in the
*peak* reached during this phase - independent of whether that freed memory
is later returned to the OS, which is exactly the point: this changes how
much is ever allocated, not how much is freed.

This composes with the already-validated `nnls_dtype='float32'` option
(`weight_solvers.py`): applying both should multiply their savings, since
they act on different axes (allocation count vs. bytes per element).

## Risk and what needs to change

- `duplicate_flip_and_interlace_orblib` and `combine_orblibs` are both
  called only from `read_vel_histograms`, per a repo-wide grep - safe to
  replace both call sites without hunting for other callers, but grep again
  before touching anything since this is a live codebase.
- Two data shapes to keep parallel: 1D histograms (`ndim==3`, losvd) and 2D
  histograms (`ndim==4`, proper motions) - the sketch above handles both,
  matching the existing branching in `duplicate_flip_and_interlace_orblib`.
- The `is_bar_disk_system()` (no-mirror) path needs its own branch in the
  merged function (sketched above as `mirror=False`) - today that case skips
  `duplicate_flip_and_interlace_orblib` entirely and calls `combine_orblibs`
  directly on the unmirrored tube data, so the merged function must
  reproduce that with a plain (non-interlaced) copy.
- `duplicate_flip_and_interlace_intmoms` (`orblib.py:1862`) is the analogous
  function for the 3D intrinsic-mass grid and has the same 2-generation
  pattern, but that array is far smaller (~0.004-0.014 GB observed vs. the
  multi-GB `vel_histograms`) - not worth touching for memory reasons, could
  be left alone or done later purely for consistency.
- The old functions should stay in the codebase temporarily during
  validation (e.g. as `_legacy_duplicate_flip_and_interlace_orblib`) so the
  new path can be checked against them array-for-array before deleting them.

## Prototype results (2026-07-31)

The `combine_and_mirror_orblibs` sketch above was implemented standalone
(not yet in `orblib.py`) and tested against the current
`duplicate_flip_and_interlace_orblib` + `combine_orblibs` pipeline on the
mirrored/1D-histogram path (`mirror=True`, `ndim==3`), the dominant case:

- **Bit-exact equivalence** (`np.array_equal`, not `allclose`) on both the
  small NGC6278 test orbit library and the 1.2GB scratch build from the OOM
  investigation.
- **Peak RSS**, measured with a background sampling thread (not
  before/after snapshots - the peak is transient and happens mid-function)
  and compared across separate fresh processes (an in-process comparison
  was tried first and gave a nonsensical result, contaminated by allocator
  history from running both paths in one process - the same effect this
  whole investigation is about):

  | orbit library | old path peak | new path peak | reduction |
  |---|---|---|---|
  | small (0.36GB final array) | 1.224 GB | 1.166 GB | 4.7% |
  | 1.2GB scratch build | 2.942 GB | 2.586 GB | 12.1% |

Not yet tested: `mirror=False` (bar/disk systems - `dev_tests/bartest.yaml`
exists but its orbit library has never been integrated) and 2D histograms
/ proper motions (`dev_tests/user_test_config_ml_with_pm.yaml` exists but
its orbit library lacks the `*_pm_hist.dat.bz2` files). Both branch
explicitly on `orblib.y.ndim`/a `mirror` flag in the sketch above and are
expected to behave the same way, but that is not yet demonstrated the way
the two tested cases are.

## Validation plan before trusting this

This is a pure refactor - same output, less peak memory - so it can be
verified by exact equivalence, not just "chi2 looks similar":

1. On the existing NGC6278 orbit library (`dev_tests/NGC6278_output`, small
   enough to iterate on quickly) and the 1.2GB scratch orbit library already
   built during the OOM investigation, run both the old path
   (`duplicate_flip_and_interlace_orblib` + `combine_orblibs`) and the new
   `combine_and_mirror_orblibs`, and assert `np.array_equal` (not just
   `allclose` - this should be bit-identical, since it's the same additions
   in the same order, just fewer intermediate copies) on the resulting
   `.y` arrays for every kinematic dataset.
2. Repeat on a bar/disk-system config if one exists in `dev_tests`, to
   exercise the `mirror=False` branch specifically.
3. Re-run the RSS instrumentation already used earlier in this
   investigation (`read_vel_histograms` phase-boundary snapshots) before and
   after, on the same 1.2GB/3.5GB/7.7GB scratch orbit libraries, to confirm
   the predicted peak reduction actually shows up - this is measurable
   independent of the OS-reclaim question, since it's about the peak
   reached, not what happens after.
4. Only after (1)-(3) pass, run a full weight-solve end to end (chi2/KKT/gap
   vs. the current code) to confirm nothing downstream changed.

## What this does and doesn't fix

Fixes: the transient 2-3x peak during the read/combine phase, for any
dataset size, without relying on the allocator returning anything.

Does not fix: the ~24% (matrix-only) to potentially ~40-50% (if applied
upstream) memory reduction available from `nnls_dtype='float32'` - these are
independent, complementary levers. Does not fix cross-model RSS creep in a
long-lived `Pool` worker (that's `ncpus_weights_maxtasksperchild`, already
implemented).
