# Potential next tasks — MaxwellBloch

A survey of the codebase after the 0.8.0 release. Tasks are grouped by effort and nature.

---

## 🐛 Quick bug fixes

| File | Issue |
|------|-------|
| ~~`bin/make-gif-ffmpeg.sh:59`~~ | ✅ Fixed |
| ~~`bin/make-mp4-fixed-frame.py:17`~~ | ✅ Fixed (`scipy.ndimage.interpolation` → `scipy.ndimage`) |
| ~~`bin/make-mp4-fixed-frame-2-fields.py:17`~~ | ✅ Fixed |
| ~~`docs/usage/scripts.md:5,27`~~ | ✅ Fixed |

---

## 📄 Documentation updates

| File | Issue |
|------|-------|
| ~~`README.md:33`~~ | ✅ Updated to `qutip>=5` |
| ~~`docs/install.md:8`~~ | ✅ Updated |
| ~~`docs/docs-environment.yml`~~ | ✅ Deleted |
| ~~`.readthedocs.yml`~~ | ✅ Updated to Python 3.11, modern config format |
| ~~`CHANGELOG.md`~~ | ✅ 0.8.0 entry added |
| ~~`docs/conf.py`~~ | ✅ Copyright year updated to 2019–2026 |

---

## 🧪 Test gaps

| Item | Detail |
|------|--------|
| ~~`ob_base.py` has no test file~~ | ✅ `test_ob_base.py` added: TestSigma, TestSetH0, TestHIList |
| ~~`TestSolveOverThermalDetunings` skipped~~ | ✅ Dead skip removed; method doesn't exist |
| ~~`test_t_funcs.py:97`~~ | ✅ `TestSech.test_fwhm` added — verifies sech FWHM matches argument |
| ~~GH#198~~ | ✅ `TestVoigtProfile.test_voigt_approaches_lorentzian_for_small_thermal_width` added |
| `test_mb_solve.py` | Many MBSolve methods untested — `TestBuildZlist`, `TestCheck` added but more needed |

---

## 🔧 Code quality

| Item | Detail |
|------|--------|
| ~~Missing docstrings~~ | ✅ `utility.py`, `ob_base.py`, `ob_solve.py` docstrings added; `t_funcs.py` factory functions annotated |
| ~~GH#106~~ | ✅ `fixed.rabi_freq_abs` removed; callers updated to `fixed.rabi_freq(..., part="abs")` |
| ~~GH#7~~ | ✅ `ramp_onoff` and `ramp_offon` already compose `ramp_on`/`ramp_off` closures |
| GH#178 | Clarify public API: use underscore convention for private methods throughout |
| ~~GH#133~~ | ✅ Vectorised `get_fields_sum_coherence` — precomputed index arrays, NumPy advanced indexing + matmul |
| ~~GH#24~~ | ✅ Split `build_velocity_classes` into `_normalise_velocity_classes` + `_build_thermal_delta_list`; vectorised `z_step_fields_euler`/`_ab` with broadcasting |
| 48 TODO/FIXME comments | Notable: `ob_base.py:117`, `mb_solve.py:598` (duplicates OBAtom), `hyperfine.py:242` (validate I,J,F range) |
| ~~Type annotations~~ | ✅ All 14 source files fully annotated: `sigma.py`, `utility.py`, `t_funcs.py`, `fixed.py`, `spectral.py`, `field.py`, `ob_base.py`, `ob_atom.py`, `ob_solve.py`, `mb_solve.py` |
| ~~Positional arguments~~ | ✅ Internal call sites updated to use keyword arguments throughout |
| Mutable default arguments | `B006` (flake8-bugbear) not in ruff `select`; many function signatures use `[]` and `{}` as defaults — enable `B` ruleset and fix violations |

---

## 🏗️ CI / infrastructure

| Item | Detail |
|------|--------|
| ~~Docs build not in CI~~ | ✅ Docs job added to `.github/workflows/ci.yml` with `nbsphinx_execute=never` |
| ~~Pre-commit lint hook~~ | ✅ `.git/hooks/pre-commit` runs ruff check + format check |
| ~~README badge wrong URL~~ | ✅ Updated from `python-package-conda.yml` to `ci.yml` |
| ~~NumPy 2 notebook compatibility~~ | ✅ Fixed `.tolist()`/`float()` repr change and `np.trapz` → `np.trapezoid` across 7 notebooks |
| ~~Benchmark regression check in CI~~ | ✅ `bench` job added to CI: runs `test_bench_two_level_sech_2pi`, caches baseline, fails on >25% mean regression |

---

## 🔬 Features / investigations

| Item | Detail |
|------|--------|
| GH#141 | QuTiP 5 has a proper `Coefficient` API — investigate whether the custom `intp()` t_func and time-dep solver wiring can be simplified |
| GH#166 | Add `macro.py` module for physical unit conversions |
| GH#157 | Add `Field.upper_levels()` / `lower_levels()` convenience methods |
| GH#151 | Add selection rule validation to `LevelF.__init__` — enforce `|J - I| <= F <= J + I` (already a TODO at `hyperfine.py:242`) |
| GH#71/#83 | Omegas refactor — 88 occurrences across 8 files; old 2017 investigation into moving Omegas methods to OBAtom |

---

## ✅ Issues closed

| Issue | Reason |
|-------|--------|
| ~~GH#7 — ramp_onoff refactor~~ | ✅ Closed — closures already compose `ramp_on`/`ramp_off` |
| ~~GH#29 — Fix PEP8 problems~~ | ✅ Closed — ruff applied across the whole codebase in 0.8.0 |
| ~~GH#106 — Remove deprecated `fixed.rabi_freq_abs`~~ | ✅ Closed — removed; callers updated to `fixed.rabi_freq(..., part="abs")` |
| ~~GH#139 — Add CHANGELOG to manifest~~ | ✅ Closed — CHANGELOG.md and MANIFEST.in both exist |
| ~~GH#163 — Write a version_info.py~~ | ✅ Closed — replaced by `importlib.metadata` + `bump-my-version` in 0.8.0 |
| ~~GH#198 — Voigt profile test~~ | ✅ Closed — `TestVoigtProfile` added in `test_spectral.py` |
| ~~GH#215 — Set up Travis deploy~~ | ✅ Closed — GitHub Actions publish workflow with OIDC in 0.8.0 |
| ~~GH#237 — Fix nbsphinx docs bug~~ | ✅ Closed — resolved in PR #240 (merged) |
