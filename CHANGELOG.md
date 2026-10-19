# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic
Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Time-grid checks in `PropagationSolver.solve()` (with `check_stability`). It warns if the fields or summed coherences
  are not resolved by `t_steps` (the error of rebuilding them from every other time sample, from a cubic spline like
  QuTiP's, is above `t_step_error_warn`; the suggested `t_steps` uses the convergence order measured between strides
  2 and 4, kept in `solver.t_step_error_estimate`), and if a field is still changing near the edge of the time window
  by more than `window_edge_warn` of its peak, so a pulse or its ringing is probably cut off
  (`solver.window_edge_estimate`). `suggest_t_steps()` estimates a grid before solving, from the input bandwidth, Rabi
  frequencies, detunings and Doppler width.
- Numerical-stability checks in `PropagationSolver.solve()`. A blow-up (non-finite field, peak field more than
  `blowup_factor` times its input, or a failed master-equation integration) now stops the run with
  `PropagationInstabilityError`, giving the z reached and a suggested `z_steps_inner`. The Euler/Adams-Bashforth
  difference gives a free estimate of the error per z-step (`solver.step_error_estimate`); `solve()` warns when it
  exceeds `step_error_warn`. `check_stability=False` turns this off. `z_steps_inner="auto"` repeats the solve with more
  inner steps until the estimate is within tolerance.
- Automatic differentiation of field propagation (optional, experimental). `PropagationSolver.to_autodiff()` returns
  a pure, `jit`-able function of a parameter pytree (pulse parameters, `rabi_freq`, interaction strengths,
  detunings, decay scale, or sampled input fields) that reproduces `solve()`, with exact gradients of any scalar loss
  of the propagated field in reverse or forward mode (`AutodiffPropagation.value_and_grad`). A gradient costs about
  five to six propagations regardless of the number of parameters. Install with `pip install 'clerq[jax]'`; JAX must run
  in float64 (`jax.config.update("jax_enable_x64", True)`). Supports several co-propagating fields, velocity
  classes, bare energies and the `gaussian`, `sech`, `square` and `ramp_*` pulses (the `on`/`off` edge times of a
  `square` pulse can be set but not differentiated, since the field is discontinuous there); counter-propagating
  fields are not supported. Importing `clerq` never imports JAX. Results are validated against `PropagationSolver`
  for both schemes and every `z_steps_inner`. Documented in the API reference and in the example notebook
  `autodiff-pulse-design`.
- `AutodiffPropagation.minimize()` runs an optax optimiser (Adam by default, or any `optax` transformation) on a loss of
  the propagated field and returns an `OptimizeResult`. `optax` is now part of the `jax` extra.

- `PropagationSolver.rabi_freqs_zt`: the solved field as a *cyclic* Rabi frequency, the same unit as `rabi_freq` in the
  problem definition. Shape `(num_fields, z_steps+1, t_steps+1)`.

### Changed
- `AutodiffPropagation.propagate()` (and so the argument of every loss function) returns the field as a cyclic Rabi
  frequency in the units of `rabi_freq` and `omega0`, equal to `PropagationSolver.rabi_freqs_zt`. Previously it was
  angular, so `propagate()` output could not be passed straight back as `omega0`. This feature is unreleased, so
  nothing outside this branch is affected.

### Deprecated
- `PropagationSolver.Omegas_zt` (the field as an *angular* frequency, which needed `/ (2 * np.pi)` to compare with
  `rabi_freq`) is deprecated in favour of `rabi_freqs_zt`. It still works and raises a `DeprecationWarning`. The
  arrays returned by `solve()`, and what `fields_area()` and the `clerq.spectral`, `clerq.fixed`, `clerq.fileio` and
  `clerq.plot` helpers compute from them, are unchanged and still angular.

### Fixed
- Calling `PropagationSolver.solve()` more than once on the same object gave different answers (for example transmitted
  energy fractions of 0.96848, 0.95962, 0.95950, ... for a sech pulse), and editing a field's `rabi_freq_t_args`
  between solves was not picked up correctly. The z-stepping replaced `system.H_Omega_list` with the last slice's
  interpolated coefficients and never restored it, so the next solve drove slice 0 with stale data instead of the
  analytic pulse. `solve()` (and the deprecated `mbsolve()`) now rebuilds the analytic `H_Omega_list` before and after
  solving, so repeated solves are identical to a fresh solver.
- `PropagationSolver.solve(step="ab")` (the default) stored every output slice one inner step early whenever
  `z_steps_inner > 1`: the Euler bootstrap recorded slice 1 after a single inner step instead of `z_steps_inner`
  steps, and every later slice inherited the offset, so `Omegas_zt[:, j]` was at
  `z_min + (j * z_steps_inner - 1) * h_inner` rather than at `zlist[j]`. The first slice now takes the Euler
  bootstrap step followed by the remaining inner steps with Adams-Bashforth, and every slice is at its nominal z.
  `z_steps_inner=1` and `step="euler"` are unchanged. The bug also read the number density one inner step ahead
  (the second inner step's density was dropped), and the two errors cancel when the first two inner steps are in
  vacuum, so a medium that starts a couple of inner steps after `z_min` (e.g. `z_min=-0.2` with the medium from
  `z=0`) gave the same results before and after. **Results change when the medium (or a non-zero density)
  reaches `z_min`**, with `z_steps_inner > 1`, which includes the default `z_steps_inner=2`; the error is first
  order in the inner step.
- The source term of the Maxwell propagation equation weighted the coherence of each coupled pair by `factor**2`
  while the Hamiltonian coupled it with `factor`, so a field coupling several transitions with non-unit factors
  responded as `sum(factor**3)` instead of the line strength `sum(factor**2)`. The weights are now the signed
  factors themselves (the dipole matrix elements), as in the Hamiltonian. Two-level problems and fields whose
  factors are all 1 are unchanged. Hyperfine and multi-transition problems change: the absorption is stronger by
  the ratio `sum(factor**2) / sum(factor**3)`, and a linearly polarised field on a `J=1 -> 2` manifold, which had
  antisymmetric factors with `sum(factor**3) = 0` and so passed through the medium unchanged, now absorbs.
  The Rb87 hyperfine example notebooks and `usage/structure` are recalibrated (interaction strengths scaled by the
  field factor, e.g. 250 -> 72 and 50 -> 14, rounded) so they show the same physical picture as before, and re-executed.
  Cached `.qu` results solved with the old weights are not reused when a field has factors other than 1: `solve()`
  recalculates them, and `load_results()` raises a `ValueError`.
- Velocity classes (Doppler broadening) detuned the Hamiltonian by 2π times too much. `thermal_delta_*` and
  `thermal_width` are cyclic frequencies like the `decay` rates and `rabi_freq`s (as documented), and the weights
  already used them that way, but the angular `thermal_delta_list` was passed to `set_H_Delta`, which expects cyclic
  detunings and multiplies by 2π again. A class at `thermal_delta = x` now detunes the field by `x`, so the Gaussian of
  width `thermal_width` is the physical lineshape. **This changes the numerical results of every problem with
  `velocity_classes`**: the old behaviour is equivalent to a `thermal_width` and `thermal_delta_*` that are 2π times
  larger, so multiply them by 2π to reproduce an old result.
- `spectral.voigt_two_linear_known` took `thermal_width` as 2π times the Gaussian width its docstring describes, which
  matched the bug above. It now agrees with a direct convolution of the Lorentzian with `exp(-(f / thermal_width)**2)`.
- Examples recalibrated: the two-pulse and stimulated photon echoes use `thermal_width = 0.94` (was `0.15`) to keep the
  same Doppler width; the other velocity-class notebooks keep their parameters, now with the Doppler widths their
  text describes.

### Changed
- The main classes have clearer names and are exported from the top-level package, so `import clerq` is enough:

  | Before | Now |
  |---|---|
  | `clerq.mb_solve.MBSolve` | `clerq.PropagationSolver` |
  | `clerq.ob_solve.OBSolve` | `clerq.LindbladSolver` |
  | `clerq.ob_atom.OBAtom` | `clerq.System` |
  | `clerq.field.Field` | `clerq.Field` |
  | `MBSolve.mbsolve()`, `OBSolve.obsolve()` | `solve()` |
  | `atom` (argument, attribute, JSON key) | `system` |

  The modules are now `clerq.propagation`, `clerq.lindblad`, `clerq.system` and `clerq.system_base` (`OBBase` is
  `SystemBase`). "System" replaces "atom" because the solvers work for any few-level open quantum system.
- Functions that take a solver name that argument `solver` everywhere: it was `mb_solve` in `spectral`, `fixed` and
  `fileio`, and `mbs` in `plot` (`mbs_list`, in `plot.spectrum_overlay`, is `solvers`).
- Problem definitions are written with a `system` key. A `PropagationSolver` `.qu` savefile made by an earlier version no longer matches
  its problem definition, so it is solved again and overwritten the first time it is used.
- `PropagationSolver` and `System` have class docstrings describing their arguments.

### Deprecated
- The old names above still work and raise a `DeprecationWarning`: the `mb_solve`, `ob_solve`, `ob_atom` and `ob_base`
  modules and their classes, `mbsolve()`, `obsolve()`, the `atom` argument and attribute, and the `atom` key in JSON
  problem definitions, and the `mb_solve`, `mbs` and `mbs_list` keyword arguments. `mbsolve_euler()` and `mbsolve_ab()` are deprecated in favour of `solve(step=...)`. They will be
  removed in a future release. The `mbsolve` and `obsolve` command-line scripts are unchanged.

### Security
- Lockfile updated to clear 19 Dependabot alerts in transitive docs/dev packages (`tornado`, `notebook`, `jupyterlab`,
  `urllib3`, `mistune`, `anyio`, `soupsieve`; `setuptools` dropped from the lock). None are dependencies of the
  installed package (`qutip`, `tqdm`), so `pip install clerq` is unaffected.

## [Unreleased]

### Fixed
- Release assets: only the wheel, sdist and their attestations are attached to the GitHub release (a stray
  `default.gitignore` from `dist/.gitignore` was included before).
- Documentation: API page titles in the browser tab keep their case (`mb_solve.MBSolve`, not `mb_solve.mbsolve`); the
  docstrings that produced 16 quartodoc warnings now match their signatures (`hyperfine`, `ob_atom`, `spectral`).

### Changed
- CHANGELOG entries for 0.10.0, 0.11.0 and 0.12.0 added.

## [0.13.1] 2026-10-09

### Fixed
- The `mbsolve` and `obsolve` command-line scripts were broken and are fixed. `obsolve` called `OBSolve.solve`, which
  was renamed `obsolve` in 0.10.0; `mbsolve` passed the `pbar_chunk_size` keyword removed in 0.11.0, and also failed
  while building its argument parser on Python 3.14 (an unescaped `%` in a help string). The dead `--pbarchunksize`
  option is removed (the progress bar is always on), and `tests/test_bin.py` now solves a small problem with each script.
- Example notebooks recalibrated after the factor-of-2 fix in 0.13.0: interaction strengths halved (and the
  `linear-absorption` labels and Beer-Lambert formula updated), so every example solves the physical problem it
  was designed for, bit for bit as in 0.12.0; plot limits widened where the data left the axes.

### Changed
- The package is now described as classical Maxwell fields propagating through open quantum systems (Maxwell's
  equations coupled to the Lindblad master equation), not as a Maxwell-Bloch solver: PyPI summary, README, docs.
- Docs: executed notebook outputs are regenerated from scratch (the site had shown tqdm progress bars on several
  pages), a rename announcement page, and a "Why a new name" section. `make docs_execute` now stores the outputs
  in the notebooks.
- CI: the benchmark regression gate is widened from 50% to 200%, because identical code varied by more than that
  between GitHub runners.

### Added
- `CITATION.cff` (software citation metadata; also read by Zenodo and GitHub's
  "Cite this repository").
- The stimulated three-pulse photon echo example.

## [0.13.0] 2026-10-09

### Fixed
- **Breaking: corrected a factor-of-2 error in the Maxwell propagation
  equation** (`1j` → `2j` in `MBSolve`), with the spectral-theory curves in
  `clerq.spectral` updated to match. Propagation previously gave half the
  physical optical depth (the Maxwell equation used the convention for coupling
  Ω, the Bloch solver the standard Ω/2), so absorption and propagation results
  change. Recompute any cached `.qu` results.

### Changed
- **Renamed from `maxwellbloch` to `clerq`.** The package, import name and PyPI
  name are now `clerq`; update `import maxwellbloch` to `import clerq`. The
  Plotly theme template is now named `"clerq"`, and `.qu` save files record the
  `clerq` version instead of the `maxwellbloch` version. Earlier entries below
  keep the old name.
- Documentation moved from Sphinx/readthedocs to Quarto (clerq.org). Executed
  executed outputs are stored in the notebooks; the `.qu` solver caches
  are no longer committed.

## [0.12.0] 2026-05-06

### Added
- Plot helpers `plot.field_profile` and `plot.field_z_profile_anim`; `show_z_bounds` and `speed_of_light` (fixed-frame
  transform) options on the plot functions; x-axis bounds for `plot.spectrum`.
- `SECH_FWHM_CONV` constant in `utility`.

### Changed
- Performance: the number-density profile N(z) is computed once up front and the sigma-operator loop is vectorised.
- `OBBase.mesolve` passes `e_ops` to `qutip.mesolve` by keyword.
- Examples consolidated: the SIT, linear-absorption, V-system and Λ EIT/slow-light notebooks are merged into one
  notebook each.

## [0.11.0] 2026-05-03

### Added
- Counter-propagating fields: `Field(counter_propagating=True)` models a field travelling against the others via a
  Doppler-shift factor. A depletion check raises `CounterPropagatingDepletionError` when forward-only propagation is no
  longer a valid approximation (disable with `check_counter_prop_depletion=False`).
- `maxwellbloch.plot`: a Plotly plotting module with a shared theme (field space-time and envelope plots, spectra,
  populations, coherences, pulse area). Install with the `plot` extra.
- `spectral` FFT functions accept a `window` (SciPy window name); `plot.spectrum` gains `freq_range`, an `arcsinh`
  frequency scale and `window`.
- Examples: Autler-Townes splitting, CPT dark resonance, adiabatons, two-pulse photon echo, odd-π sech pulses,
  Rb 87 D2 cycling transition, and a Rydberg EIT draft (kept out of the published docs); an examples gallery page.
- `bench_max` benchmark (Rb 87 V-simulton, D1 + D2 hyperfine, thermal).

### Changed
- **Breaking:** the print-based progress reporting of `MBSolve.mbsolve` and its `pbar_chunk_size` argument are
  replaced by a tqdm progress bar, controlled by `progress` (default `True`) and `progress_show_area`.
- `spectral.absorption`, `dispersion` and `rabi_freq` take the optional `window` argument.

## [0.10.0] 2026-04-25

### Added
- Savefile integrity check (GH#199): `.qu` files store a hash of the problem definition and the library versions;
  loading a file built from a different problem raises `ValueError`, and older files load with a warning.
- Generated API reference pages; a `docs_serve` Makefile target.

### Changed
- **Breaking:** `OBSolve.solve` is renamed `OBSolve.obsolve`.
- `maxwell_boltzmann` moved to `utility` and gained an `offset` parameter.
- Hyperfine classes share their JSON (de)serialisation through a mixin (`to_json_str`, `get_json_dict`); dead code
  removed and in-code TODOs resolved.

### Fixed
- A missing savefile directory is created before the `.qu` file is written (GH#195).
- `IndexError` in `get_fields_sum_coherence` when fields are added after construction.

## [0.9.0] 2026-04-24

### Added
- `Field.upper_levels()` and `Field.lower_levels()` — return sorted unique upper/lower level indices for a field (GH#157)
- `LevelF` now validates that `|J−I| ≤ F ≤ J+I`, raising `ValueError` for invalid hyperfine quantum numbers (GH#151)

### Changed
- **Performance**: replaced `intp()` closure (which recreated a `scipy.interpolate.interp1d` on every ODE function call) with QuTiP 5 `InterCoefficient` built once per z-step — ~43× speedup on the inner propagation loop (GH#141)
- Internal methods across `ob_atom.py`, `ob_solve.py`, `mb_solve.py`, and `field.py` renamed with `_` prefix to reflect their private status (GH#178); e.g. `build_H_0` → `_build_H_0`, `build_c_ops` → `_build_c_ops`
- `_z_step_fields_euler` and `_z_step_fields_ab` now accept pre-computed `h` and `N` scalars instead of z-coordinates; `z_prev` (unused) removed from AB signature (GH#71)
- `coherences()` in `ob_atom.py` vectorised with precomputed index arrays; stale per-field loop removed
- `build_velocity_classes` refactored; `_z_step_fields` vectorised across fields (GH#24)
- Pre-allocate `_Omegas_z_buf` in `init_Omegas_zt`; inner loop reuses buffer via `np.copyto` (GH#71)

### Tests
- Added unit tests for `_normalise_velocity_classes` and `_build_thermal_delta_list`
- Added `TestBuildIntpHOmegaList` for `_build_intp_H_Omega_list`
- Added `TestZStepFields` for `_z_step_fields_euler` and `_z_step_fields_ab`
- Added tests for `Field.upper_levels()`, `Field.lower_levels()`, and `LevelF` validation

## [0.8.1] 2026-04-23

- Added type annotations across all source files
- Fixed mutable default arguments (replaced `[]`/`{}` defaults with `None`)
- Enabled flake8-bugbear (`B`) ruleset in ruff; fixed all B006/B007 violations
- Added docs build job to CI; fixed NumPy 2 compatibility in notebooks
- Added `test_ob_base.py`; added tests for sech FWHM and Voigt profile
- Removed deprecated `fixed.rabi_freq_abs`; updated callers to `fixed.rabi_freq(..., part="abs")`
- Used keyword arguments at all internal call sites
- Updated `.readthedocs.yml` to Python 3.11 and modern config format

## [0.8.0] 2026-04-15

- Updated to support QuTiP 5, NumPy 2, and SciPy 1.14+
- Migrated packaging from `setup.py` to `pyproject.toml` (PEP 517/518/621)
- Switched primary environment manager from conda to uv; added `uv.lock`
- Added ruff for linting and formatting
- Replaced Travis CI with GitHub Actions; added OIDC trusted publishing to PyPI
- Added `bump-my-version` for semver release management
- Dropped Python < 3.10 support

## [0.7.1] 2020-08-10

- Fixed bug on adding multiple fields

## [0.7.0] 2020-07-11

- Added `mbsolve` and `obsolve` scripts
- Added ability to set sech pulse width by defining the full-width at half max
- Added ability to set sech and gaussian pulse amplitudes by defining the area
    as multiples of pi
- Changed the default Rabi frequency of a field from 0.0 to 1.0
- Added docs on Structure and Angular Momentum
- Added Examples docs

## [0.6.0] 2020-05-17

- Fixed Doppler broadening
- Added docs on Usage, Examples

## [0.5.0] 2019-08-07

- Added hyperfine level strength factor method
- Added hyperfine factors method for isotropic field (equal components in all 
    three polarisations)
- Added methods for summing populations and coherences coupled by a field
- Fixed the detuning term of the interaction Hamiltonian for fields coupling
    multiple upper levels
- Added functionality to build structure for single valence electron atoms
- Added ability to set QuTiP options
- Added factors to fields and atom decays for atomic structure
- Added ability to set initial state of atom
- Added solver options to OBSolve
- Added decay factors representing structure in coupled levels
- Added factors list to Field object for strength factors on coupled levels
- Added function to calculate hyperfine Clebsch-Gordan coefficients
- Added Wigner 3j and Wigner 6j functions
- Added script to make MP4s of systems with two fields

## [0.4.0] 2018-08-28

- Fixed the mp4 and gif filenames so they don't need to be .json.mp4.gif
- Changed t_funcs so that the _1, _2, … suffix is no longer required
- Changed ob_atom to atom
- Added label to atom
- Fixed interpolation on mp4 script
- Changed look of animated plots
- Fixed mp4 script to use real part of rabi freq
- Added automatic update of version number from setup.py
- Changed to allow empty velocity class to be set
- Changed the time function for fields to be constant if none is given

## [0.3.0] 2018-01-08

- Added Travis CI build
- Fixed bug that tests submodule not in package
- Added git hash to unreleased version number
- Added scripts to make MP4s and gifs to /bin
- Added methods to save field CSV files
- Refactored spectral methods to separate module
- Refactored fixed frame methods to separate module
- Fixed bug where empty decay list would cause exception
- Fixed bug where fixed frame didn't work if no inner z steps
- Added testing.py to run tests
- Added link to video and gif tools
- Fixed multiple coupled levels bug
- Added MBSolve.set_field_rabi_freq_t_func() and set_field_rabi_freq_t_args()
    methods so we can add custom input fields
- Added OBAtom.build() method for resetting operators
- Added FFT methods for spectral analysis of pulse propagation
- Added ability to run `ob_solve` and `mb_solve` from the command line
- Added methods to get results in the fixed frame of reference
- Changed to QuTiP progress bar to remove a dependency

## [0.2.0] 2017-03-27

- Added ability to solve the Maxwell-Bloch equations
