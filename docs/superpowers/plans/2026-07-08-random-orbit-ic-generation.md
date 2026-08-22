# Random (Agama-backed) Orbit IC Generation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in orbit initial-condition (IC) generator, backed by Agama's Jeans-equation sampling, as a drop-in replacement for the Fortran `orbitstart`/`orbitstart_bar` binaries — zero changes to the Fortran orbit integrator, classification, mirroring, or weight-solving.

**Architecture:** A new `dynamite/random_ic.py` module builds an `agama.Potential` from a model's existing MGE/dark-halo/BH parameters, samples tube-family (full 3D velocity) and box-family (`v_y=0`) seed points via Agama's density sampling + Jeans-equation velocity ellipsoid, dithers each seed into `dithering**3` realizations, and writes `begin.dat`/`beginbox.dat` in the exact fixed-width format the Fortran binary already produces. A new `orblib_settings.ic_generator` config flag (`grid` default, `random` opt-in) selects this path inside `OrbitLibrary.get_orbit_ics()` in `dynamite/orblib.py`.

**Tech Stack:** Python 3, NumPy, Agama (`pip install agama`), pytest. Reference fixture: a real `begin.dat` already on disk at `dev_tests/NGC6278_bayesopt_qml_output/models/orblib_002_001/datfil/begin.dat`.

## Global Constraints

- Zero changes to any `legacy_fortran/*.f90` file.
- Zero changes to orbit classification (`integrator_find_orbtype`), mirroring, or weight-solving code paths.
- Default config behavior (`ic_generator: grid`, or the setting absent) must be byte-for-byte unchanged from today.
- The `noreg` field is written as a constant `0` for every orbit — confirmed dead code downstream, not reproduced.
- Validation target is `dev_tests/user_test_config_ml.yaml` (nE=6, nI2=5, nI3=4, dithering=1 → 120 orbits), not a full ω Cen run.
- Box-family orbits get `v_y = 0` forced (matching the existing a priori box-seeding convention); everything else is sampled from Agama's Jeans ellipsoid.

---

## File Structure

```
dynamite/random_ic.py                          new — potential translation, sampling, dithering, file writing
dynamite/config_reader.py                       modify — add ic_generator setting + validation
dynamite/orblib.py                              modify — branch in get_orbit_ics()
dev_tests/test_random_ic/
    conftest.py                                 shared fixtures (small test config, real begin.dat fixture)
    fixtures/begin_dat_reference.dat             copy of a real orbitstart-produced begin.dat
    test_ic_format.py                            byte-for-byte round-trip test against the reference fixture
    test_potential_translation.py                Agama potential vs dynamite's own potential evaluation
    test_sampling.py                             position/velocity sampling shape + physical sanity checks
    test_ic_writer.py                             dithering + file-writing tests
    test_end_to_end.py                            small-config grid vs random comparison
dev_tests/user_test_config_ml_random_ic.yaml     new — copy of user_test_config_ml.yaml with ic_generator: random
```

**Interface contract:**
- **Consumes:** the same physical parameters `OrbitLibrary.create_fortran_input_orblib()` already has access to on `self` — `self.stars` (MGE + kinematic/population data), `self.system.get_all_dark_non_plummer_components()` (dark halo), `self.system.get_component_from_class(physys.Plummer)` (BH), `self.parset`, `self.settings` (orblib_settings dict: `nE`, `nI2`, `nI3`, `dithering`, `logrmin`, `logrmax`, `random_seed`).
- **Produces:** `{mod_dir}/datfil/begin.dat`, `{mod_dir}/datfil/beginbox.dat`, in the same format `orbitstart`/`orbitstart_bar` write today, consumed unmodified by `OrbitLibrary.get_orbit_library()`.

---

## Task 1: Confirm Agama installs and its basic API works

**Files:**
- Create: `dev_tests/test_random_ic/test_agama_smoke.py`

No dynamite integration yet — this task only confirms Agama is usable in the target dev environment (the `gc-diag` conda env at `/opt/miniconda3/envs/gc-diag`, per this project's existing environment notes) before any other task depends on it.

- [ ] **Step 1: Install Agama in the gc-diag environment**

```bash
conda run -n gc-diag pip install agama
```

If this fails, STOP and report the failure — every later task in this plan depends on Agama being installable here. Do not proceed to Task 2 until this succeeds.

- [ ] **Step 2: Write the smoke test**

```python
# dev_tests/test_random_ic/test_agama_smoke.py
import numpy as np


def test_agama_imports():
    import agama
    agama.setUnits(length=1, velocity=1, mass=1)  # arbitrary internal units for this smoke test


def test_agama_builds_plummer_potential_and_evaluates_density():
    import agama
    agama.setUnits(length=1, velocity=1, mass=1)
    pot = agama.Potential(type='Plummer', mass=1.0, scaleRadius=1.0)
    # density at the center of a unit-mass, unit-scale-radius Plummer sphere
    # is a known finite positive value
    rho0 = pot.density(0.0, 0.0, 0.0)
    assert rho0 > 0
    assert np.isfinite(rho0)


def test_agama_circular_velocity_positive():
    import agama
    agama.setUnits(length=1, velocity=1, mass=1)
    pot = agama.Potential(type='Plummer', mass=1.0, scaleRadius=1.0)
    vcirc = pot.circularVelocity(1.0)
    assert vcirc > 0
    assert np.isfinite(vcirc)


def test_agama_sample_returns_points_and_masses():
    import agama
    agama.setUnits(length=1, velocity=1, mass=1)
    pot = agama.Potential(type='Plummer', mass=1.0, scaleRadius=1.0)
    xyz, mass = pot.sample(100)
    assert xyz.shape == (100, 3)
    assert mass.shape == (100,)
    assert np.all(np.isfinite(xyz))
```

- [ ] **Step 3: Run the smoke tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_agama_smoke.py -v
```
Expected: all 4 tests PASS. If `pot.sample` doesn't exist under this name in the installed Agama version, check `import agama; help(agama.Potential)` and record the actual sampling API name — this affects Task 5's implementation, so resolve it here before moving on.

- [ ] **Step 4: Commit**

```bash
git add dev_tests/test_random_ic/test_agama_smoke.py
git commit -m "test: confirm Agama installs and basic API works in gc-diag env"
```

---

## Task 2: Lock down the exact IC file format against a real fixture

**Files:**
- Create: `dev_tests/test_random_ic/fixtures/begin_dat_reference.dat`
- Create: `dev_tests/test_random_ic/test_ic_format.py`
- Create: `dynamite/random_ic.py` (format constants only, for now)

**Interfaces:**
- Produces: `dynamite.random_ic.IC_LINE_FORMAT` (a format-string constant), `dynamite.random_ic.parse_ic_line(line: str) -> dict` (parses one fixed-width line into a dict with keys `iE, iI2, iI3, x, y, z, vx, vy, vz, rcirc, tcirc, vcirc, noreg`), `dynamite.random_ic.format_ic_line(iE, iI2, iI3, x, y, z, vx, vy, vz, rcirc, tcirc, vcirc, noreg) -> str`.

A real `begin.dat` already exists on disk at `dev_tests/NGC6278_bayesopt_qml_output/models/orblib_002_001/datfil/begin.dat` — copy it in as a reference fixture rather than generating a fresh one.

- [ ] **Step 1: Copy the reference fixture**

```bash
mkdir -p dev_tests/test_random_ic/fixtures
cp dev_tests/NGC6278_bayesopt_qml_output/models/orblib_002_001/datfil/begin.dat \
   dev_tests/test_random_ic/fixtures/begin_dat_reference.dat
```

- [ ] **Step 2: Write the failing test**

```python
# dev_tests/test_random_ic/test_ic_format.py
import os

FIXTURE = os.path.join(os.path.dirname(__file__), 'fixtures', 'begin_dat_reference.dat')


def test_parse_header_line():
    from dynamite.random_ic import parse_ic_header
    with open(FIXTURE) as f:
        header = f.readline()
    nE, nI2, nI3 = parse_ic_header(header)
    assert nE * nI2 * nI3 == 24  # 24 data lines follow in this fixture
    assert (nE, nI2, nI3) == (2, 4, 3)


def test_parse_first_data_line_matches_known_values():
    from dynamite.random_ic import parse_ic_line
    with open(FIXTURE) as f:
        f.readline()  # skip header
        first_line = f.readline()
    parsed = parse_ic_line(first_line)
    assert parsed['iE'] == 1
    assert parsed['iI2'] == 1
    assert parsed['iI3'] == 1
    assert abs(parsed['x'] - 6.1906350331e+14) / 6.1906350331e+14 < 1e-10
    assert parsed['y'] == 0.0
    assert abs(parsed['z'] - 3.1122423983e+15) / 3.1122423983e+15 < 1e-10
    assert parsed['noreg'] == 0


def test_format_then_parse_round_trip():
    from dynamite.random_ic import parse_ic_line, format_ic_line
    with open(FIXTURE) as f:
        f.readline()
        original_line = f.readline()
    parsed = parse_ic_line(original_line)
    rewritten = format_ic_line(**parsed)
    reparsed = parse_ic_line(rewritten)
    assert parsed == reparsed


def test_full_fixture_round_trip_byte_identical():
    """Parsing every data line and reformatting it must reproduce the
    original file byte-for-byte, proving format compatibility with the
    real orbitstart output."""
    from dynamite.random_ic import parse_ic_line, format_ic_line
    with open(FIXTURE) as f:
        lines = f.readlines()
    header, data_lines = lines[0], lines[1:]
    rewritten_lines = [format_ic_line(**parse_ic_line(l)) for l in data_lines]
    for original, rewritten in zip(data_lines, rewritten_lines):
        assert original.rstrip('\n') == rewritten.rstrip('\n'), (
            f"mismatch:\n  original:  {original!r}\n  rewritten: {rewritten!r}"
        )
```

- [ ] **Step 2: Run test to verify it fails**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_ic_format.py -v
```
Expected: `ImportError` — `dynamite.random_ic` does not exist yet.

- [ ] **Step 3: Implement the format constants and parse/format functions**

The fixed-width format is `(3I5, 9ES30.10, I4)`: three 5-character integers, nine 30-character scientific-notation floats, one 4-character integer. Fortran `ES` format uses a leading `0.` mantissa style with capital or lowercase `E`; the reference fixture uses uppercase `E` with a `+`/`-` sign and 10-digit mantissa, e.g. `6.1906350331E+14`.

```python
# dynamite/random_ic.py
"""Agama-backed orbit initial-condition generation.

Drop-in replacement for the Fortran orbitstart/orbitstart_bar binaries:
produces begin.dat/beginbox.dat in the exact same fixed-width format,
so the unmodified Fortran orbit integrator (orblib/orblibbox) can read
them without any changes.
"""

# Column widths matching Fortran format (3I5, 9ES30.10, I4)
_INT_WIDTH = 5
_FLOAT_WIDTH = 30
_FLOAT_PRECISION = 10
_NOREG_WIDTH = 4

_FLOAT_FIELDS = ('x', 'y', 'z', 'vx', 'vy', 'vz', 'rcirc', 'tcirc', 'vcirc')
_INT_FIELDS = ('iE', 'iI2', 'iI3')


def parse_ic_header(header_line):
    """Parse the 'nE nI2 nI3' header line of an IC file.

    Returns
    -------
    (nE, nI2, nI3) : tuple of int
    """
    parts = header_line.split()
    return int(parts[0]), int(parts[1]), int(parts[2])


def _format_fortran_float(value):
    """Format a float as Fortran ES30.10: sign, mantissa, E, exponent sign+2 digits."""
    formatted = f"{value:.10E}"
    mantissa, exponent = formatted.split('E')
    exp_sign = '+' if exponent[0] != '-' else '-'
    exp_digits = exponent.lstrip('+-').zfill(2)
    if not mantissa.startswith('-'):
        mantissa = mantissa  # keep positive sign implicit, matches fixture (no leading '+')
    return f"{mantissa}E{exp_sign}{exp_digits}"


def parse_ic_line(line):
    """Parse one fixed-width IC data line into a dict of field values."""
    pos = 0
    result = {}
    for name in _INT_FIELDS:
        result[name] = int(line[pos:pos + _INT_WIDTH])
        pos += _INT_WIDTH
    for name in _FLOAT_FIELDS:
        result[name] = float(line[pos:pos + _FLOAT_WIDTH])
        pos += _FLOAT_WIDTH
    result['noreg'] = int(line[pos:pos + _NOREG_WIDTH])
    return result


def format_ic_line(iE, iI2, iI3, x, y, z, vx, vy, vz, rcirc, tcirc, vcirc, noreg):
    """Format one IC data line matching Fortran's (3I5, 9ES30.10, I4)."""
    parts = [f"{iE:5d}", f"{iI2:5d}", f"{iI3:5d}"]
    for value in (x, y, z, vx, vy, vz, rcirc, tcirc, vcirc):
        parts.append(f"{_format_fortran_float(value):>{_FLOAT_WIDTH}}")
    parts.append(f"{noreg:{_NOREG_WIDTH}d}")
    return ''.join(parts)
```

- [ ] **Step 4: Run tests, fix formatting mismatches**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_ic_format.py -v
```
Expected: all 4 tests PASS. If `test_full_fixture_round_trip_byte_identical` fails on whitespace/sign details, inspect the mismatch message (it prints both strings) and adjust `_format_fortran_float`/`format_ic_line` until it's byte-identical — do not relax the test.

- [ ] **Step 5: Commit**

```bash
git add dev_tests/test_random_ic/fixtures/begin_dat_reference.dat \
        dev_tests/test_random_ic/test_ic_format.py \
        dynamite/random_ic.py
git commit -m "feat: lock IC file fixed-width format against real begin.dat fixture"
```

---

## Task 3: Add the `ic_generator` config setting

**Files:**
- Modify: `dynamite/config_reader.py:103-109` (the `for quad in ['nr', 'nth', 'nph']` defaulting block in `Settings.validate()`)
- Test: `dev_tests/test_random_ic/test_config_setting.py`

**Interfaces:**
- Produces: `config.settings.orblib_settings['ic_generator']`, one of `'grid'` or `'random'`, defaulting to `'grid'` if absent from the YAML.

- [ ] **Step 1: Write the failing test**

```python
# dev_tests/test_random_ic/test_config_setting.py
import pytest


def test_ic_generator_defaults_to_grid(tmp_path):
    from dynamite.config_reader import Settings
    import logging
    settings = Settings(
        logger=logging.getLogger('test'),
        orblib_settings={
            'nE': 6, 'logrmin': -1.0, 'logrmax': 2.0, 'nI2': 5, 'nI3': 4,
            'dithering': 1,
        },
        parameter_space_settings={'generator_type': 'GridWalk',
                                  'which_chi2': 'kinchi2'},
        io_settings={},
        weight_solver_settings={'type': 'NNLS', 'nnls_solver': 'scipy'},
        multiprocessing_settings={},
    )
    settings.validate()
    assert settings.orblib_settings['ic_generator'] == 'grid'


def test_ic_generator_rejects_invalid_value():
    from dynamite.config_reader import Settings
    import logging
    settings = Settings(
        logger=logging.getLogger('test'),
        orblib_settings={
            'nE': 6, 'logrmin': -1.0, 'logrmax': 2.0, 'nI2': 5, 'nI3': 4,
            'dithering': 1, 'ic_generator': 'nonsense',
        },
        parameter_space_settings={'generator_type': 'GridWalk',
                                  'which_chi2': 'kinchi2'},
        io_settings={},
        weight_solver_settings={'type': 'NNLS', 'nnls_solver': 'scipy'},
        multiprocessing_settings={},
    )
    with pytest.raises(ValueError):
        settings.validate()
```

- [ ] **Step 2: Run test to verify it fails**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_config_setting.py -v
```
Expected: FAIL — `settings.orblib_settings['ic_generator']` raises `KeyError` (setting doesn't exist yet), and the invalid-value test fails because nothing raises.

- [ ] **Step 3: Add the setting to `Settings.validate()`**

In `dynamite/config_reader.py`, immediately after the existing `quad_nr`/`quad_nth`/`quad_nph` defaulting block (lines 103-109):

```python
        for quad in ['nr', 'nth', 'nph']:
            key = 'quad_' + quad
            if key not in self.orblib_settings.keys():
                default = 10 if quad == 'nr' else 6
                self.orblib_settings[key] = default
                self.logger.info(f'No value given for orblib setting {key} '
                                 f'- set to its default {default}.')
        if 'ic_generator' not in self.orblib_settings:
            self.orblib_settings['ic_generator'] = 'grid'
            self.logger.info("No value given for orblib setting "
                             "ic_generator - set to its default 'grid'.")
        if self.orblib_settings['ic_generator'] not in ('grid', 'random'):
            text = "orblib_settings: ic_generator must be 'grid' or " \
                   f"'random', but is {self.orblib_settings['ic_generator']!r}."
            self.logger.error(text)
            raise ValueError(text)
```

- [ ] **Step 4: Run tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_config_setting.py -v
```
Expected: both tests PASS.

- [ ] **Step 5: Commit**

```bash
git add dynamite/config_reader.py dev_tests/test_random_ic/test_config_setting.py
git commit -m "feat: add ic_generator config setting (grid|random, default grid)"
```

---

## Task 4: Agama potential translation

**Verified against Agama's own reference implementation** (`py/schwarzlib.py`
and `py/example_forstand.py` in github.com/GalacticDynamics-Oxford/Agama —
the actual Forstand code) rather than guessed API names. Key corrections
this made to the original draft of this task:
- There is no `type='MGE'` Agama potential. Each MGE Gaussian component is
  built as a separate `agama.Density(density='Spheroid', ...)` with a
  specific mass/axis-ratio/cutoff recipe (`schwarzlib.getDensityParamsMGE`),
  and all components are summed via `agama.Density(*components)`
  (`schwarzlib.makeDensityMGE`).
- The combined stellar+halo density is wrapped into one potential solver
  via `agama.Potential(type='Multipole', density=combined_density, lmax=.., mmax=.., gridSizeR=..)`,
  not a per-component potential sum.
- The BH is added separately as a `Plummer` potential and combined via
  `agama.Potential(pot_gal, pot_bh)` (composing two `Potential` objects,
  since Multipole expansion isn't needed/suitable for a point mass).

**Files:**
- Modify: `dynamite/random_ic.py`
- Test: `dev_tests/test_random_ic/test_potential_translation.py`

**Interfaces:**
- Consumes: `stars.mge_pot.data` (structured array with fields `I`, `sigma`, `q`, `PA_twist` — confirmed in `dynamite/mges.py:69`), a dark-halo component object (from `system.get_all_dark_non_plummer_components()`), a BH component object (from `system.get_component_from_class(physys.Plummer)`), viewing angles `(theta, phi, psi)` in radians, `self.parset['ml']`.
- Produces: `dynamite.random_ic.mge_to_agama_density(mge_pot_data, distance_kpc, inclination_rad) -> agama.Density`, `dynamite.random_ic.build_agama_potential(mge_pot_data, dh_density, bh_mass, bh_scale, theta, phi, psi, ml, distMPc) -> agama.Potential`.

This task only needs to produce *a* valid combined potential; validating it matches dynamite's own Fortran potential to high precision is out of scope for this plan (that level of cross-validation belongs to the end-to-end task, where chi2 comparison is the real test). Here we validate structural correctness: the returned object evaluates finite, positive density and circular velocity at a few representative radii.

- [ ] **Step 1: Write the failing test**

```python
# dev_tests/test_random_ic/test_potential_translation.py
import numpy as np


def _fake_mge_pot_data():
    # 2-Gaussian MGE, columns I (Lsun/pc^2), sigma (arcsec), q, PA_twist (deg)
    dtype = [('I', 'f8'), ('sigma', 'f8'), ('q', 'f8'), ('PA_twist', 'f8')]
    return np.array([(1000.0, 1.0, 0.8, 0.0),
                     (200.0, 5.0, 0.7, 0.0)], dtype=dtype)


def test_mge_to_agama_density_total_mass_is_positive_and_finite():
    from dynamite.random_ic import mge_to_agama_density
    density = mge_to_agama_density(_fake_mge_pot_data(), distance_kpc=10000.0,
                                   inclination_rad=np.pi / 2)
    mass = density.totalMass()
    assert np.isfinite(mass) and mass > 0


def test_build_agama_potential_stellar_only_is_finite_and_positive():
    from dynamite.random_ic import build_agama_potential, mge_to_agama_density
    stellar_density = mge_to_agama_density(_fake_mge_pot_data(),
                                           distance_kpc=10000.0,
                                           inclination_rad=np.pi / 2)
    pot = build_agama_potential(
        stellar_density=stellar_density,
        dh_density=None,
        bh_mass=0.0,
        bh_scale=0.0,
    )
    rho = pot.density(1.0, 0.0, 0.0)
    vcirc = pot.circularVelocity(1.0)
    assert np.isfinite(rho) and rho > 0
    assert np.isfinite(vcirc) and vcirc > 0


def test_build_agama_potential_with_bh_increases_inner_vcirc():
    from dynamite.random_ic import build_agama_potential, mge_to_agama_density
    stellar_density = mge_to_agama_density(_fake_mge_pot_data(),
                                           distance_kpc=10000.0,
                                           inclination_rad=np.pi / 2)
    pot_no_bh = build_agama_potential(
        stellar_density=stellar_density, dh_density=None,
        bh_mass=0.0, bh_scale=0.0,
    )
    pot_with_bh = build_agama_potential(
        stellar_density=stellar_density, dh_density=None,
        bh_mass=1e8, bh_scale=1e-4,
    )
    # Adding a BH should raise circular velocity very close to the centre
    r_inner = 1e-3
    assert pot_with_bh.circularVelocity(r_inner) > pot_no_bh.circularVelocity(r_inner)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_potential_translation.py -v
```
Expected: `ImportError` — `mge_to_agama_density`/`build_agama_potential` don't exist yet.

- [ ] **Step 3: Implement `mge_to_agama_density` and `build_agama_potential`**

This follows `schwarzlib.getDensityParamsMGE`/`makeDensityMGE` exactly
(verified against github.com/GalacticDynamics-Oxford/Agama's
`py/schwarzlib.py`), adapted to dynamite's MGE column names (`I`, `sigma`,
`q` — dynamite's `mges.py:69` — instead of Agama's example script's raw
`tab[:,0]`/`tab[:,1]`/`tab[:,2]` array columns, but the same underlying
recipe).

```python
# append to dynamite/random_ic.py
import numpy as np
import agama


def _mge_gaussian_to_spheroid_params(mass, Sx, Sy, Sz):
    """One MGE Gaussian component's Agama Density parameters.

    Matches schwarzlib.getDensityParamsMGE: approximates a 3D Gaussian
    via Agama's generic 'Spheroid' density profile with gamma=beta=0,
    alpha=1 (flat core, no outer power-law tail) and a Gaussian-matched
    cutoff (outerCutoffRadius = sqrt(2)*Sx, cutoffStrength=2).
    """
    return dict(
        density='Spheroid',
        axisRatioY=Sy / Sx,
        axisRatioZ=Sz / Sx,
        scaleRadius=1,
        gamma=0,
        beta=0,
        alpha=1,
        outerCutoffRadius=2 ** 0.5 * Sx,
        cutoffStrength=2,
        densityNorm=mass / ((2 * np.pi) ** 1.5 * Sx * Sy * Sz),
    )


def mge_to_agama_density(mge_pot_data, distance_kpc, inclination_rad):
    """Deproject dynamite's observed MGE into an Agama Density object.

    Parameters
    ----------
    mge_pot_data : structured array with fields 'I' (Lsun/pc^2), 'sigma'
        (arcsec), 'q' (observed axis ratio) — dynamite's ``stars.mge_pot.data``
    distance_kpc : float — distance to the system [kpc]
    inclination_rad : float — inclination angle (dynamite's viewing angle
        theta) [radians]

    Returns
    -------
    agama.Density
    """
    arcsec_to_kpc = distance_kpc * np.pi / 648000.0
    if 1 - np.min(mge_pot_data['q']) ** 2 > np.sin(inclination_rad) ** 2:
        raise ValueError('Deprojection is impossible for the given '
                         'inclination — check theta/q consistency.')
    components = []
    for row in mge_pot_data:
        I_row, sigma_arcsec, q_row = row['I'], row['sigma'], row['q']
        Sx = sigma_arcsec * arcsec_to_kpc
        Sy = Sx
        Sz = Sx * (1 - (1 - q_row ** 2) / np.sin(inclination_rad) ** 2) ** 0.5
        mass = 2 * np.pi * I_row * (1000 * arcsec_to_kpc) ** 2 * q_row
        components.append(_mge_gaussian_to_spheroid_params(mass, Sx, Sy, Sz))
    return agama.Density(*components)


def build_agama_potential(stellar_density, dh_density, bh_mass, bh_scale):
    """Compose the full Agama Potential: stars + optional dark halo,
    solved via a Multipole expansion, plus an optional BH added directly.

    Parameters
    ----------
    stellar_density : agama.Density (from mge_to_agama_density)
    dh_density : agama.Density or None (from Task 4b's NFW translation)
    bh_mass : float — BH mass [Msun], 0 to omit
    bh_scale : float — BH softening/scale radius, should be tiny (e.g. 1e-4
        in the system's length unit) relative to the smallest orbit radii
        of interest — not a physical size, just numerical softening

    Returns
    -------
    agama.Potential
    """
    combined_density = agama.Density(stellar_density, dh_density) \
        if dh_density is not None else stellar_density
    pot_galaxy = agama.Potential(type='Multipole', density=combined_density,
                                 lmax=8, mmax=0, gridSizeR=25)
    if bh_mass and bh_mass > 0:
        pot_bh = agama.Potential(type='Plummer', mass=bh_mass,
                                 scaleRadius=max(bh_scale, 1e-8))
        return agama.Potential(pot_galaxy, pot_bh)
    return pot_galaxy
```

> `lmax=8, mmax=0, gridSizeR=25` are much smaller than Forstand's example
> defaults (`lmax=32, mmax=0 or 6, gridSizeR=40`), matching this plan's
> small validation target (`dev_tests/user_test_config_ml.yaml`'s tiny
> 120-orbit library) rather than a production-scale galaxy model. `mmax=0`
> assumes an axisymmetric or spherical system; revisit if the validation
> config's system is triaxial (`mmax=6`, matching Forstand's triaxial
> case, per `example_forstand.py`'s `mmax=0 if symmetry[0]!='t' else 6`).

- [ ] **Step 4: Run tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_potential_translation.py -v
```
Expected: all PASS. If any Agama parameter name above doesn't match the
installed version (Agama's API has had minor renames across versions),
check `help(agama.Density)`/`help(agama.Potential)` in the gc-diag env
and adjust — record what you find as a comment in the code.

- [ ] **Step 5: Commit**

```bash
git add dynamite/random_ic.py dev_tests/test_random_ic/test_potential_translation.py
git commit -m "feat: build Agama potential from stellar MGE + optional BH"
```

---

## Task 4b: Dark halo translation (NFW only, matching the validation config)

**Verified against Agama's own reference implementation**
(`schwarzlib.makeDensityNFWHalo`). Important correction from the original
draft of this task: Agama's literal `type='NFW'` potential has **infinite
total mass** and cannot be used directly for density sampling (`.sample()`
needs a finite-mass density to draw from) — Forstand instead builds NFW as
a truncated `'spheroid'` **Density** (`alpha=1, beta=3, gamma=1`, i.e. the
classic NFW power-law slopes, cut off at large radius so the mass
converges), not a `Potential` object, and it gets summed into the same
Multipole expansion as the stellar density in `build_agama_potential`
(Task 4).

The validation target (`dev_tests/user_test_config_ml.yaml`) needs to be
checked for which dark-halo type it actually uses before this task
starts — do not assume NFW.

- [ ] **Step 1: Check which dark-halo type the validation config uses, and its parameter names**

```bash
grep -B2 -A8 -i "dh\|halo\|nfw" dev_tests/user_test_config_ml.yaml
```
Record the component's `type` field and parameter names. If it's `NFW`,
proceed with Steps 2-4 below. Then check `dynamite/physical_system.py`'s
`NFW` class for its `par_names` attribute and, critically, its Fortran
acceleration formula (`legacy_fortran/dmpotent.f90`'s `dm_accel` NFW
case, already read earlier this session: uses `rho_c`/`r_c` with
`t1 = -4*pi*G*rho_c*r_c**3/d**2`) — dynamite's NFW is parametrized by
central density `rho_c` and scale radius `r_c`, **not** by peak circular
velocity like Agama's `makeDensityNFWHalo(rscale, vcirc)` expects, so a
conversion is needed (Step 3 below).

- [ ] **Step 2: Write the failing test**

```python
# append to dev_tests/test_random_ic/test_potential_translation.py
def test_nfw_translation_has_finite_mass_and_positive_vcirc():
    from dynamite.random_ic import nfw_to_agama_density
    # rho_c [Msun/pc^3], r_c [pc] — arbitrary physically-reasonable values
    density = nfw_to_agama_density(rho_c=0.01, r_c=1000.0)
    mass = density.totalMass()
    assert np.isfinite(mass) and mass > 0


def test_nfw_density_positive_at_scale_radius():
    from dynamite.random_ic import nfw_to_agama_density
    density = nfw_to_agama_density(rho_c=0.01, r_c=1000.0)
    rho_at_rc = density.density(1000.0, 0.0, 0.0)
    assert np.isfinite(rho_at_rc) and rho_at_rc > 0
```

- [ ] **Step 2: Run test to verify it fails**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_potential_translation.py::test_nfw_translation_has_finite_mass_and_positive_vcirc -v
```
Expected: `ImportError` — `nfw_to_agama_density` doesn't exist yet.

- [ ] **Step 3: Implement `nfw_to_agama_density`**

Converts dynamite's `(rho_c, r_c)` parametrization directly into Agama's
truncated-spheroid NFW density norm, without going through Agama's
`vcirc`-based helper (avoids an unnecessary round-trip conversion, and
matches dynamite's own Fortran normalization directly): dynamite's NFW
density profile is `rho(r) = rho_c / [(r/r_c)(1+r/r_c)^2]`, which is
exactly the density implied by `alpha=1, beta=3, gamma=1` in Agama's
generic Spheroid profile with `densityNorm=rho_c` and `scaleRadius=r_c`
— no conversion needed for these two parameters, only the truncation
radius must be chosen (default: `100 * r_c`, matching Forstand's
`makeDensityNFWHalo` default, since the outer truncation only affects
convergence far outside the region orbits are sampled from).

```python
# append to dynamite/random_ic.py

def nfw_to_agama_density(rho_c, r_c, rcutoff=None):
    """Build an Agama Density for an NFW halo from dynamite's own
    (rho_c, r_c) parametrization (legacy_fortran/dmpotent.f90's NFW case).

    Parameters
    ----------
    rho_c : float — NFW characteristic density [Msun/pc^3]
    r_c : float — NFW scale radius [pc]
    rcutoff : float or None — outer truncation radius (Agama's literal NFW
        potential has infinite mass, so this Spheroid is truncated for
        finite-mass density sampling); default 100*r_c, matching Forstand.

    Returns
    -------
    agama.Density
    """
    if rcutoff is None:
        rcutoff = 100 * r_c
    return agama.Density(type='spheroid', alpha=1, beta=3, gamma=1,
                         densitynorm=rho_c, scaleradius=r_c,
                         outercutoffradius=rcutoff)
```

Then in the code that assembles the full potential (Task 7's
`generate_ic_files`, not `build_agama_potential` itself — `dh_density` is
already a parameter of `build_agama_potential` from Task 4), call this
function to produce that `dh_density` argument:

```python
    dh_density = None
    if dh_component is not None:
        dh_density = nfw_to_agama_density(
            rho_c=dh_parset_values['rho_c'], r_c=dh_parset_values['r_c'],
        )
```

(Adjust the parset key names `rho_c`/`r_c` to match whatever Step 1 found
in `NFW.par_names` — do not guess; use the actual attribute names.)

- [ ] **Step 4: Run tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_potential_translation.py -v
```
Expected: all PASS. If Agama's `Spheroid` density parameter names
(`densitynorm`/`scaleradius`/`outercutoffradius`, lowercase, matching
`schwarzlib.py`'s exact usage) don't match the installed version, check
`help(agama.Density)` in the gc-diag env and adjust.

- [ ] **Step 5: Commit**

```bash
git add dynamite/random_ic.py dev_tests/test_random_ic/test_potential_translation.py
git commit -m "feat: translate dynamite NFW dark halo to truncated Agama Spheroid density"
```

---

## Task 5: Position and velocity sampling

**Verified against Agama's own reference implementation**
(`example_forstand.py:576-580`). Important correction from the original
draft of this task: Agama's `Density.sample(N, potential=..., beta=...,
kappa=...)` does position **and** velocity sampling in one native call —
`beta` is the velocity anisotropy (0 = isotropic, up to ~0.5), `kappa` is
rotation sense (`1`/`-1`) — solving the axisymmetric Jeans equations
internally. There is no need to hand-roll a Jeans solver or an isotropic
approximation; the hand-rolled version from the original draft of this
task is dropped entirely in favor of this real API.

**Files:**
- Modify: `dynamite/random_ic.py`
- Test: `dev_tests/test_random_ic/test_sampling.py`

**Interfaces:**
- Consumes: `agama.Density` (the same combined stellar+halo density built
  in Task 4/4b, *not* the composite `Potential` — Agama's `.sample()` is a
  `Density` method, called with the full `Potential` passed as its
  `potential=` argument so it can solve the Jeans equations self-
  consistently against gravity from all components), the corresponding
  `agama.Potential`, `n_bundles: int`, `box_fraction: float` (fraction of
  bundles seeded as box-family via post-hoc `vy=0`, default `0.2` — an
  open tuning parameter per the design doc, not derived from Agama), an
  anisotropy `beta: float` (default `0.0`, isotropic — dynamite's targets
  like ω Cen are not disky/rotating, so `kappa`/rotation is not used
  here, unlike Forstand's disk-galaxy example), `random_seed: int`.
- Produces: `dynamite.random_ic.sample_seed_points(density, potential, n_bundles, box_fraction, beta, random_seed) -> dict` with keys `x, y, z, vx, vy, vz` (each shape `(n_bundles,)`) and `is_box` (boolean array, shape `(n_bundles,)`).

- [ ] **Step 1: Write the failing test**

```python
# dev_tests/test_random_ic/test_sampling.py
import numpy as np
import agama


def _test_density_and_potential():
    agama.setUnits(length=1, velocity=1, mass=1)
    density = agama.Density(type='Plummer', mass=1.0, scaleRadius=1.0)
    potential = agama.Potential(type='Multipole', density=density,
                                lmax=0, mmax=0, gridSizeR=20)
    return density, potential


def test_sample_seed_points_shapes():
    from dynamite.random_ic import sample_seed_points
    density, potential = _test_density_and_potential()
    seeds = sample_seed_points(density, potential, n_bundles=50,
                               box_fraction=0.2, beta=0.0, random_seed=42)
    for key in ('x', 'y', 'z', 'vx', 'vy', 'vz', 'is_box'):
        assert seeds[key].shape == (50,)
    assert np.all(np.isfinite(seeds['x']))
    assert np.all(np.isfinite(seeds['vx']))


def test_box_family_has_zero_vy():
    from dynamite.random_ic import sample_seed_points
    density, potential = _test_density_and_potential()
    seeds = sample_seed_points(density, potential, n_bundles=50,
                               box_fraction=0.2, beta=0.0, random_seed=42)
    assert np.all(seeds['vy'][seeds['is_box']] == 0.0)
    assert np.any(seeds['vy'][~seeds['is_box']] != 0.0)


def test_box_fraction_is_approximately_respected():
    from dynamite.random_ic import sample_seed_points
    density, potential = _test_density_and_potential()
    seeds = sample_seed_points(density, potential, n_bundles=1000,
                               box_fraction=0.2, beta=0.0, random_seed=42)
    observed_fraction = np.mean(seeds['is_box'])
    assert abs(observed_fraction - 0.2) < 0.05


def test_reproducible_with_same_seed():
    from dynamite.random_ic import sample_seed_points
    density, potential = _test_density_and_potential()
    seeds1 = sample_seed_points(density, potential, n_bundles=50,
                                box_fraction=0.2, beta=0.0, random_seed=7)
    seeds2 = sample_seed_points(density, potential, n_bundles=50,
                                box_fraction=0.2, beta=0.0, random_seed=7)
    np.testing.assert_array_equal(seeds1['x'], seeds2['x'])
```

- [ ] **Step 2: Run test to verify it fails**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_sampling.py -v
```
Expected: `ImportError` — `sample_seed_points` doesn't exist yet.

- [ ] **Step 3: Implement `sample_seed_points`**

```python
# append to dynamite/random_ic.py

def sample_seed_points(density, potential, n_bundles, box_fraction, beta,
                       random_seed):
    """Sample orbit seed positions and velocities via Agama's native
    axisymmetric-Jeans density sampling.

    A `box_fraction` of seeds get vy forced to zero after sampling
    (box-family), matching dynamite's existing a priori box-orbit
    convention — Agama/Forstand itself has no such split (confirmed:
    no "box" orbit family anywhere in Agama's schwarzlib.py or
    example_forstand.py), this is purely to keep feeding dynamite's
    existing two-Fortran-integrator-binary architecture unchanged.

    Parameters
    ----------
    density : agama.Density — the combined stellar+halo density to draw
        positions from (same object passed into build_agama_potential)
    potential : agama.Potential — full potential (stars+halo+BH), used
        by Agama internally to solve the Jeans equations self-
        consistently
    n_bundles : int — number of distinct orbit bundles to seed
    box_fraction : float in [0, 1) — fraction seeded as box-family
    beta : float — velocity anisotropy passed to Agama's Jeans solver
        (0 = isotropic)
    random_seed : int

    Returns
    -------
    dict with keys x, y, z, vx, vy, vz, is_box (each shape (n_bundles,))
    """
    agama.setRandomSeed(random_seed)  # Agama has its own RNG (see
                                      # example_forstand.py's agama.setRandomSeed)
    points, _ = density.sample(n_bundles, potential=potential, beta=beta)
    x, y, z, vx, vy, vz = points.T

    rng = np.random.default_rng(random_seed)
    is_box = rng.random(n_bundles) < box_fraction
    vy = np.where(is_box, 0.0, vy)

    return {'x': x, 'y': y, 'z': z, 'vx': vx, 'vy': vy, 'vz': vz,
           'is_box': is_box}
```

- [ ] **Step 4: Run tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_sampling.py -v
```
Expected: all PASS. If `density.sample(..., potential=..., beta=...)`
returns a different column order or shape than `(N, 6)` in the installed
Agama version, check `help(agama.Density.sample)` in the gc-diag env —
`example_forstand.py` calls it as `densityStars.sample(N, potential=pot_fidu, beta=0.3, kappa=1)[0]`
and treats the result as directly usable as `ic` (an Nx6 orbit initial-
conditions array) for `agama.schwarzlib.runModel`, so `(N, 6)` with
columns `(x, y, z, vx, vy, vz)` is the expected shape — adjust the
unpacking above only if this doesn't hold in practice.

- [ ] **Step 5: Commit**

```bash
git add dynamite/random_ic.py dev_tests/test_random_ic/test_sampling.py
git commit -m "feat: sample orbit seed positions/velocities via Agama's native Jeans sampling"
```

---

## Task 6: Dithering, circular-orbit quantities, and IC file writing

**Files:**
- Modify: `dynamite/random_ic.py`
- Test: `dev_tests/test_random_ic/test_ic_writer.py`

**Interfaces:**
- Consumes: output of `sample_seed_points` (Task 5), `agama.Potential` (for `circularVelocity`), `dithering: int`.
- Produces: `dynamite.random_ic.dither_seeds(seeds, dithering, random_seed) -> dict` (same keys as `sample_seed_points`, but each array expanded from shape `(n_bundles,)` to `(n_bundles * dithering**3,)`), `dynamite.random_ic.write_ic_file(path, seeds, potential, nE, nI2, nI3)` (writes one IC file).

- [ ] **Step 1: Write the failing test**

```python
# dev_tests/test_random_ic/test_ic_writer.py
import numpy as np
import agama
import tempfile
import os


def _test_potential():
    agama.setUnits(length=1, velocity=1, mass=1)
    return agama.Potential(type='Plummer', mass=1.0, scaleRadius=1.0)


def test_dither_seeds_expands_shape():
    from dynamite.random_ic import dither_seeds
    seeds = {
        'x': np.array([1.0, 2.0]), 'y': np.array([0.0, 0.0]),
        'z': np.array([0.0, 0.0]), 'vx': np.array([0.1, 0.2]),
        'vy': np.array([0.1, 0.0]), 'vz': np.array([0.0, 0.0]),
        'is_box': np.array([False, True]),
    }
    dithered = dither_seeds(seeds, dithering=2, random_seed=1)
    for key in ('x', 'y', 'z', 'vx', 'vy', 'vz', 'is_box'):
        assert dithered[key].shape == (2 * 8,)  # 2 seeds * 2**3 dithers


def test_dither_preserves_box_vy_zero():
    from dynamite.random_ic import dither_seeds
    seeds = {
        'x': np.array([1.0]), 'y': np.array([0.0]), 'z': np.array([0.0]),
        'vx': np.array([0.1]), 'vy': np.array([0.0]), 'vz': np.array([0.0]),
        'is_box': np.array([True]),
    }
    dithered = dither_seeds(seeds, dithering=3, random_seed=1)
    assert np.all(dithered['vy'] == 0.0)


def test_write_ic_file_produces_correct_line_count_and_header():
    from dynamite.random_ic import dither_seeds, write_ic_file
    pot = _test_potential()
    seeds = {
        'x': np.array([1.0, 2.0]), 'y': np.array([0.0, 0.0]),
        'z': np.array([0.0, 0.0]), 'vx': np.array([0.1, 0.2]),
        'vy': np.array([0.1, 0.0]), 'vz': np.array([0.0, 0.0]),
        'is_box': np.array([False, False]),
    }
    dithered = dither_seeds(seeds, dithering=2, random_seed=1)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'begin.dat')
        write_ic_file(path, dithered, pot, nE=2, nI2=2, nI3=2)
        with open(path) as f:
            lines = f.readlines()
    header = lines[0].split()
    assert int(header[0]) * int(header[1]) * int(header[2]) == 16  # 2*8
    assert len(lines) == 1 + 16


def test_write_ic_file_lines_are_parseable():
    from dynamite.random_ic import (dither_seeds, write_ic_file,
                                    parse_ic_line)
    pot = _test_potential()
    seeds = {
        'x': np.array([1.0]), 'y': np.array([0.0]), 'z': np.array([0.0]),
        'vx': np.array([0.1]), 'vy': np.array([0.1]), 'vz': np.array([0.0]),
        'is_box': np.array([False]),
    }
    dithered = dither_seeds(seeds, dithering=1, random_seed=1)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'begin.dat')
        write_ic_file(path, dithered, pot, nE=1, nI2=1, nI3=1)
        with open(path) as f:
            f.readline()
            data_line = f.readline()
    parsed = parse_ic_line(data_line)
    assert abs(parsed['x'] - 1.0) / 1.0 < 1e-9
    assert parsed['noreg'] == 0
```

- [ ] **Step 2: Run test to verify it fails**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_ic_writer.py -v
```
Expected: `ImportError` — `dither_seeds`/`write_ic_file` don't exist yet.

- [ ] **Step 3: Implement `dither_seeds` and `write_ic_file`**

```python
# append to dynamite/random_ic.py

_DITHER_POS_SCALE = 0.1   # kpc-equivalent internal-unit spread; adjust to match rcirc scale
_DITHER_VEL_FRAC = 0.05   # fraction of local circular velocity


def dither_seeds(seeds, dithering, random_seed):
    """Expand each seed into dithering**3 nearby realisations.

    Perturbs position by a small uniform offset and velocity by a small
    fraction of the seed's own velocity magnitude — matching the "small
    uniform grid centred on its initial 6D phase-space coordinates"
    approach SchwarMAX/Forstand use, rather than dynamite's classic
    heavier per-grid-cell dithering.

    Parameters
    ----------
    seeds : dict with keys x, y, z, vx, vy, vz, is_box (each shape (n,))
    dithering : int
    random_seed : int

    Returns
    -------
    dict, same keys, each shape (n * dithering**3,)
    """
    rng = np.random.default_rng(random_seed)
    n = len(seeds['x'])
    n_dith = dithering ** 3

    result = {key: np.repeat(seeds[key], n_dith) for key in seeds}

    pos_offsets = rng.uniform(-_DITHER_POS_SCALE, _DITHER_POS_SCALE,
                              size=(n * n_dith, 3))
    result['x'] = result['x'] + pos_offsets[:, 0]
    result['y'] = result['y'] + pos_offsets[:, 1]
    result['z'] = result['z'] + pos_offsets[:, 2]

    vel_mag = np.repeat(
        np.sqrt(seeds['vx']**2 + seeds['vy']**2 + seeds['vz']**2), n_dith,
    )
    vel_mag = np.where(vel_mag > 0, vel_mag, 1.0)  # avoid zero-scale noise
    vel_offsets = rng.normal(0.0, _DITHER_VEL_FRAC, size=(n * n_dith, 3))
    result['vx'] = result['vx'] + vel_offsets[:, 0] * vel_mag
    result['vy'] = np.where(result['is_box'], 0.0,
                           result['vy'] + vel_offsets[:, 1] * vel_mag)
    result['vz'] = result['vz'] + vel_offsets[:, 2] * vel_mag

    return result


def write_ic_file(path, seeds, potential, nE, nI2, nI3):
    """Write one IC file (begin.dat or beginbox.dat) in Fortran fixed-width
    format, sequentially grouped so the Fortran integrator's own bundling
    (every dithering**3 consecutive lines = one bundle) lines up correctly.

    Parameters
    ----------
    path : str — output file path
    seeds : dict with keys x, y, z, vx, vy, vz (post-dithering, each
        shape (n_total,), already in the sequential bundle order
        produced by dither_seeds)
    potential : agama.Potential — used to compute rcirc/vcirc/tcirc
    nE, nI2, nI3 : int — header integers; their product must equal
        len(seeds['x'])
    """
    n_total = len(seeds['x'])
    assert nE * nI2 * nI3 == n_total, (
        f"header product {nE}*{nI2}*{nI3}={nE*nI2*nI3} != n_total={n_total}"
    )

    r = np.sqrt(seeds['x']**2 + seeds['y']**2 + seeds['z']**2)
    vcirc = np.array([potential.circularVelocity(ri) for ri in r])
    rcirc = r
    tcirc = 2.0 * np.pi * rcirc / np.where(vcirc > 0, vcirc, 1.0)

    with open(path, 'w') as f:
        f.write(f"{nE:12d}{nI2:12d}{nI3:12d}\n")
        for i in range(n_total):
            iE = i // (nI2 * nI3) + 1
            iI2 = (i // nI3) % nI2 + 1
            iI3 = i % nI3 + 1
            line = format_ic_line(
                iE=iE, iI2=iI2, iI3=iI3,
                x=seeds['x'][i], y=seeds['y'][i], z=seeds['z'][i],
                vx=seeds['vx'][i], vy=seeds['vy'][i], vz=seeds['vz'][i],
                rcirc=rcirc[i], tcirc=tcirc[i], vcirc=vcirc[i], noreg=0,
            )
            f.write(line + '\n')
```

- [ ] **Step 4: Run tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_ic_writer.py -v
```
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add dynamite/random_ic.py dev_tests/test_random_ic/test_ic_writer.py
git commit -m "feat: dithering and fixed-width IC file writing"
```

---

## Task 7: Wire into `OrbitLibrary.get_orbit_ics()`

**Files:**
- Modify: `dynamite/orblib.py:470-498`

**Interfaces:**
- Consumes: everything from Tasks 4-6, plus `self.settings['ic_generator']`, `self.settings['nE']`, `self.settings['nI2']`, `self.settings['nI3']`, `self.settings['dithering']`, `self.settings['random_seed']`, `self.stars`, `self.system`, `self.parset`.

- [ ] **Step 1: Write the failing test**

```python
# dev_tests/test_random_ic/test_orblib_wiring.py
from unittest.mock import MagicMock, patch


def test_get_orbit_ics_calls_random_ic_when_configured():
    """When ic_generator == 'random', get_orbit_ics must call
    dynamite.random_ic.generate_ic_files instead of shelling out to the
    orbitstart Fortran binary."""
    from dynamite import orblib

    fake_self = MagicMock()
    fake_self.settings = {'ic_generator': 'random', 'nE': 6, 'nI2': 5,
                          'nI3': 4, 'dithering': 1, 'random_seed': 42}
    fake_self.mod_dir = '/tmp/fake_model/'

    with patch('dynamite.random_ic.generate_ic_files') as mock_generate:
        with patch('subprocess.run') as mock_subprocess:
            orblib.OrbitLibrary.get_orbit_ics(fake_self)
            mock_generate.assert_called_once()
            mock_subprocess.assert_not_called()


def test_get_orbit_ics_calls_fortran_when_grid():
    """When ic_generator == 'grid' (the default), get_orbit_ics must
    keep calling the existing Fortran subprocess path unchanged."""
    from dynamite import orblib
    import os

    fake_self = MagicMock()
    fake_self.settings = {'ic_generator': 'grid'}
    fake_self.mod_dir = '/tmp/fake_model/'
    fake_self.legacy_directory = '/fake/legacy/dir'
    fake_self.system.is_bar_disk_system.return_value = False
    fake_p = MagicMock()
    fake_p.stdout.decode.return_value = ''
    fake_p.returncode = 0

    with patch('dynamite.random_ic.generate_ic_files') as mock_generate:
        with patch('subprocess.run', return_value=fake_p) as mock_subprocess:
            with patch('os.chdir'):
                orblib.OrbitLibrary.get_orbit_ics(fake_self)
                mock_subprocess.assert_called_once()
                mock_generate.assert_not_called()
```

- [ ] **Step 2: Run test to verify it fails**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_orblib_wiring.py -v
```
Expected: FAIL — `get_orbit_ics` unconditionally calls the Fortran subprocess today, so `mock_generate.assert_called_once()` fails in the first test.

- [ ] **Step 3: Branch in `get_orbit_ics`**

In `dynamite/orblib.py`, replace the `get_orbit_ics` method (currently lines 470-498):

```python
    def get_orbit_ics(self):
        """Calculate orbit initial conditions.

        Uses the Fortran orbitstart/orbitstart_bar binary by default, or
        dynamite.random_ic's Agama-backed sampler if orblib_settings'
        ic_generator is set to 'random'.
        """
        if self.settings.get('ic_generator', 'grid') == 'random':
            from dynamite import random_ic
            self.logger.info(f'Calculating initial conditions (random/Agama)'
                             f' for {self.mod_dir}.')
            random_ic.generate_ic_files(self)
            return
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        self.logger.info(f'Calculating initial conditions for {self.mod_dir}.')
        cmd = self.legacy_directory
        bar = '_bar' if self.system.is_bar_disk_system() else ''
        cmd += f'/orbitstart{bar} < infil/orbstart.in >> datfil/orbstart.log'
        p = subprocess.run(cmd,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        os.chdir(cur_dir)
        log_file = f'Logfile: {self.mod_dir}datfil/orbstart.log.'
        if not p.stdout.decode("UTF-8"):
            self.logger.info(f'...done - orbitstart{bar} exit code '
                             f'{p.returncode}. {log_file}')
        else:
            text = f'...failed! orbitstart{bar} exit code {p.returncode}. ' \
                   f'Message: {p.stdout.decode("UTF-8")}. {log_file}'
            if p.returncode == 127: # command not found
                text += 'Check DYNAMITE legacy_fortran executables.'
                self.logger.error(text)
                raise FileNotFoundError(text)
            else:
                text += f'{log_file} Be wary: DYNAMITE may crash...'
                self.logger.warning(text)
                raise RuntimeError(text)
```

Then add `generate_ic_files` to `dynamite/random_ic.py` as the single
entry point that ties Tasks 4-6 together, given an `OrbitLibrary`
instance:

```python
# append to dynamite/random_ic.py

def generate_ic_files(orblib_instance):
    """Generate begin.dat/beginbox.dat using the Agama-backed sampler.

    Parameters
    ----------
    orblib_instance : dynamite.orblib.OrbitLibrary
        Must already have .stars, .system, .parset, .settings, .mod_dir
        set (true for any OrbitLibrary at the point get_orbit_ics is
        called).
    """
    from dynamite import physical_system as physys
    stars = orblib_instance.stars
    system = orblib_instance.system
    parset = orblib_instance.parset
    settings = orblib_instance.settings

    bh = system.get_component_from_class(physys.Plummer)
    dh_list = system.get_all_dark_non_plummer_components()
    dh_component = dh_list[0] if dh_list else None

    q = parset[f'q-{stars.name}']
    p = parset[f'p-{stars.name}']
    u = parset[f'u-{stars.name}']
    theta, psi, phi = stars.triax_pqu2tpp(p, q, u)

    stellar_density = mge_to_agama_density(
        stars.mge_pot.data, distance_kpc=system.distMPc * 1e3,
        inclination_rad=theta,
    )
    dh_density = None
    if dh_component is not None:
        dh_density = nfw_to_agama_density(
            rho_c=parset[f'rho_c-{dh_component.name}'],
            r_c=parset[f'r_c-{dh_component.name}'],
        )  # adjust parset key names to match Task 4b's Step 1 findings
    combined_density = agama.Density(stellar_density, dh_density) \
        if dh_density is not None else stellar_density
    potential = build_agama_potential(
        stellar_density=stellar_density,
        dh_density=dh_density,
        bh_mass=parset[f"m-{bh.name}"],
        bh_scale=parset[f"a-{bh.name}"],
    )

    n_bundles = (settings['nE'] * settings['nI2'] * settings['nI3']) \
        // (settings['dithering'] ** 3)
    box_fraction = settings.get('ic_generator_box_fraction', 0.2)
    random_seed = settings['random_seed']

    seeds = sample_seed_points(combined_density, potential, n_bundles,
                               box_fraction, beta=0.0,
                               random_seed=random_seed)
    box_only = {k: v[seeds['is_box']] for k, v in seeds.items()}
    tube_only = {k: v[~seeds['is_box']] for k, v in seeds.items()}

    dithered_tube = dither_seeds(tube_only, settings['dithering'],
                                random_seed)
    dithered_box = dither_seeds(box_only, settings['dithering'],
                               random_seed + 1)

    n_tube = len(dithered_tube['x'])
    n_box = len(dithered_box['x'])
    write_ic_file(orblib_instance.mod_dir + 'datfil/begin.dat',
                 dithered_tube, potential, nE=n_tube, nI2=1, nI3=1)
    write_ic_file(orblib_instance.mod_dir + 'datfil/beginbox.dat',
                 dithered_box, potential, nE=n_box, nI2=1, nI3=1)
```

- [ ] **Step 4: Run tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_orblib_wiring.py -v
```
Expected: both PASS.

- [ ] **Step 5: Run the full existing dynamite test suite to check nothing else broke**

```bash
conda run -n gc-diag python -m pytest dev_tests/ -k "not test_random_ic" -v --timeout=600
```
Expected: same pass/fail state as before this task (this task must not change behavior for `ic_generator: grid`, which is every existing test's config).

- [ ] **Step 6: Commit**

```bash
git add dynamite/orblib.py dynamite/random_ic.py dev_tests/test_random_ic/test_orblib_wiring.py
git commit -m "feat: wire random_ic into OrbitLibrary.get_orbit_ics()"
```

---

## Task 8: End-to-end validation on the small dev_tests config

**Files:**
- Create: `dev_tests/user_test_config_ml_random_ic.yaml` (copy of `dev_tests/user_test_config_ml.yaml` with `ic_generator: random` added to `orblib_settings`)
- Create: `dev_tests/test_random_ic/test_end_to_end.py`

- [ ] **Step 1: Create the random-IC test config**

```bash
cp dev_tests/user_test_config_ml.yaml dev_tests/user_test_config_ml_random_ic.yaml
```

Edit `dev_tests/user_test_config_ml_random_ic.yaml`'s `orblib_settings`
block to add `ic_generator: random` (keep `nE: 6`, `nI2: 5`, `nI3: 4`,
`dithering: 1` unchanged, so both configs produce the same total
orbit count — 120 — for a fair comparison).

- [ ] **Step 2: Write the end-to-end comparison test**

```python
# dev_tests/test_random_ic/test_end_to_end.py
import os
import numpy as np
import pytest

CONFIG_DIR = os.path.dirname(os.path.dirname(__file__))


@pytest.mark.slow
def test_random_ic_produces_valid_orbit_library(tmp_path):
    """Runs the small dev_tests config with ic_generator: random end to
    end (orbit ICs -> Fortran integration -> orbit classification) and
    checks the resulting library isn't degenerate."""
    from dynamite import config_reader

    config_path = os.path.join(CONFIG_DIR, 'user_test_config_ml_random_ic.yaml')
    c = config_reader.Configuration(config_path, reset_existing_output=True)
    model = c.all_models.get_model_from_row(0)
    orblib = model.get_orblib()

    orblib.read_orbit_property_file()
    orbtypes = orblib.orb_properties  # includes classification per orbit

    # Not everything should collapse to "chaotic" (type 5) — a healthy
    # library has a real spread of X/Y/Z-tube and box classifications.
    type_counts = np.bincount(orbtypes['orbtype'].astype(int).ravel())
    n_chaotic = type_counts[5] if len(type_counts) > 5 else 0
    n_total = orbtypes['orbtype'].size
    assert n_chaotic < 0.5 * n_total, (
        f"{n_chaotic}/{n_total} orbits classified chaotic — sampling may "
        "be producing unstable/degenerate orbits"
    )


@pytest.mark.slow
def test_random_ic_chi2_is_comparable_to_grid(tmp_path):
    """Runs both the grid and random configs end to end and compares
    final chi2 — random sampling should not be dramatically worse."""
    from dynamite import config_reader

    grid_config_path = os.path.join(CONFIG_DIR, 'user_test_config_ml.yaml')
    random_config_path = os.path.join(CONFIG_DIR,
                                      'user_test_config_ml_random_ic.yaml')

    c_grid = config_reader.Configuration(grid_config_path,
                                         reset_existing_output=True)
    c_random = config_reader.Configuration(random_config_path,
                                           reset_existing_output=True)

    model_grid = c_grid.all_models.get_model_from_row(0)
    model_random = c_random.all_models.get_model_from_row(0)

    orblib_grid = model_grid.get_orblib()
    orblib_random = model_random.get_orblib()

    weight_solver_grid = model_grid.get_weights(orblib_grid)
    weight_solver_random = model_random.get_weights(orblib_random)

    chi2_grid = model_grid.chi2
    chi2_random = model_random.chi2

    # This is a proof-of-concept sanity check, not a strict quality bar:
    # random sampling should be in the same order of magnitude, not off
    # by orders of magnitude, which would indicate a real bug rather
    # than expected sampling-noise differences.
    assert chi2_random < 10 * chi2_grid, (
        f"random-IC chi2 ({chi2_random}) is more than 10x grid chi2 "
        f"({chi2_grid}) — investigate before treating this as viable"
    )
```

- [ ] **Step 3: Run the end-to-end tests**

```bash
conda run -n gc-diag python -m pytest dev_tests/test_random_ic/test_end_to_end.py -v -m slow --timeout=1200
```
Expected: both tests PASS. If `test_random_ic_produces_valid_orbit_library`
fails on excessive chaotic-orbit fraction, revisit Task 5's velocity
sampling (the isotropic `sigma ~ vcirc/sqrt(3)` approximation may be too
crude) before assuming the whole approach doesn't work. If chi2 is wildly
worse, check the Agama potential translation (Task 4/4b) first — a
mismatched potential between IC generation and the (unchanged) Fortran
integrator's own potential would explain a bad fit regardless of sampling
quality.

- [ ] **Step 4: Commit**

```bash
git add dev_tests/user_test_config_ml_random_ic.yaml dev_tests/test_random_ic/test_end_to_end.py
git commit -m "test: end-to-end validation of random IC generation on small config"
```

---

## Plan self-review notes

- **Spec coverage:** Task 1 covers the design's "confirm Agama installs" risk. Task 2 covers "re-verify IC format byte-for-byte." Tasks 4/4b/5/6 cover the full data-flow section of the design (potential build, position/velocity sampling incl. box `v_y=0`, dithering, `rcirc`/`vcirc`/`tcirc`, file writing). Task 3 covers the config flag. Task 7 covers the architecture's "invoked at the exact point orbitstart is called today" requirement. Task 8 covers the design's testing/validation section on the named dev_tests config. `noreg`-is-dead-code is handled by simply never trying to reproduce it (constant 0), consistent with the design.
- **API verification pass:** Tasks 4, 4b, and 5 were rewritten after checking Agama's actual open-source reference implementation (github.com/GalacticDynamics-Oxford/Agama, `py/schwarzlib.py` and `py/example_forstand.py` — the real Forstand code) rather than relying on guessed API names. This corrected three real mistakes in the first draft: there is no `type='MGE'` or usable `type='NFW'` potential in Agama (both need to go through `Density` objects built from `'Spheroid'` profile components instead, per the exact recipes in `getDensityParamsMGE`/`makeDensityNFWHalo`), and — most importantly — Agama's `Density.sample(N, potential=..., beta=..., kappa=...)` already performs the full axisymmetric Jeans-equation velocity draw natively, so Task 5 no longer needs (and no longer contains) a hand-rolled isotropic-Gaussian approximation. The plan now uses the real, tested API throughout, rather than leaving "verify against docs" placeholders for the biggest physics-adjacent piece of the feature.
- **Remaining open item, not a gap in this plan:** the exact conversion between dynamite's `NFW` parset keys and the `rho_c`/`r_c` names assumed in Task 4b/7's code is explicitly deferred to Task 4b's Step 1 (checking the actual `NFW.par_names`/config against the validation YAML) rather than guessed — if the real keys differ, only that one substitution needs to change, not the surrounding structure.
