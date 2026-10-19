# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

clerq is a Python package for numerically solving the propagation of classical electromagnetic fields through media of open quantum systems: Maxwell's equations coupled to the Lindblad master equation. It is used for near-resonant light in thermal atomic vapors, and models two, three, and many-level systems using density matrix formalism via the QuTiP library. (The package was formerly called MaxwellBloch.)

## Commands

```bash
# Create/sync environment (first time or after dependency changes)
uv sync --extra dev

# Run all tests in parallel
uv run pytest -n auto

# Run a single test file
uv run pytest tests/test_propagation.py

# Run a single test
uv run pytest tests/test_propagation.py::TestClassName::test_method_name

# Run with coverage
uv run pytest --cov -n auto

# Run performance benchmarks
uv run pytest tests/bench_propagation.py --benchmark-only

# Lint
uv run ruff check .

# Format
uv run ruff format .

# Build docs
make docs_html        # API pages + quarto render (uses the outputs stored in the notebooks)
make docs_execute     # re-run every notebook and store fresh outputs in the .ipynb files (~8 min)

# Build distribution
uv build

# Release — bump version, commit, tag, then push to trigger CI publish
uv run bump-my-version bump patch   # bug fixes
uv run bump-my-version bump minor   # new features, backwards-compatible
uv run bump-my-version bump major   # breaking changes
git push && git push --tags
```

Make targets: `test`, `test_cov`, `bench`, `lint`, `format`, `format_check`, `docs_html`, `docs_execute`, `docs_serve`, `dist`, `bump_patch`, `bump_minor`, `bump_major`, `clean_qu` (removes cached `.qu` files).

## Architecture

### Class hierarchy

```
SystemBase (system_base.py)              — QuTiP interface, operators, density matrix math
  └── System (system.py)                 — quantum system: energy levels, decay channels, coupled fields

LindbladSolver (lindblad.py)             — time-domain master equation solver (wraps qutip.mesolve)
  └── PropagationSolver (propagation.py) — spatial propagation through a medium

Field (field.py)                         — a laser field coupling two levels (Rabi freq, detuning)
```

`System`, `Field`, `LindbladSolver` and `PropagationSolver` are exported from the top-level `clerq` package.

LindbladSolver holds a `System` instance as `self.system` and delegates Hamiltonian construction to it.

Hyperfine structure lives in `hyperfine.py` and `angmom.py` (3j/6j symbols, Clebsch-Gordan coefficients).

### Deprecated names

Up to 0.13.1 the classes were `OBBase`, `OBAtom`, `OBSolve` and `MBSolve` (modules `ob_base`, `ob_atom`, `ob_solve`, `mb_solve`), the solve methods were `obsolve()` and `mbsolve()`, and the system was called `atom` (constructor argument, attribute and JSON key). Functions taking a solver named the argument `mb_solve` (`spectral`, `fixed`, `fileio`) or `mbs` (`plot`); it is now `solver`. These all still work and raise a `DeprecationWarning`; the helpers are in `_deprecation.py` and the tests in `tests/test_deprecated.py`. Do not use them in new code. The `mbsolve` and `obsolve` command-line scripts keep their names.

`PropagationSolver.Omegas_zt` (the field as an angular frequency) is also deprecated, in favour of the cyclic `rabi_freqs_zt`; the angular array is now the private `_Omegas_zt`, and the `solve()` return value and helpers such as `fields_area`, `spectral`, `fixed`, `fileio` and `plot` remain angular. Its warning is tested in `tests/test_deprecated.py`.

### Solver flow

- **LindbladSolver** calls `qutip.mesolve` to evolve the system density matrix in time for a single spatial point (`solve()`).
- **PropagationSolver** wraps LindbladSolver and propagates the field through space using finite differences (Adams or Euler method). At each z-step its `solve()`: (1) solves the master equation, (2) reads off coherences, (3) updates the Rabi frequency for the next step.
- Velocity classes are used for Doppler broadening: the master equation is solved independently per velocity class and results are summed.

### Configuration

Problems are defined as JSON (see `fileio.py`). `PropagationSolver` and `LindbladSolver` accept a JSON string or a dict. Test fixtures live in `tests/json/`.

### Caching

Solved results are pickled by QuTiP into `.qu` files next to the script. Use `make clean_qu` to remove them.

### Time functions

`t_funcs.py` provides named callables (sech, Gaussian, square) for pulse shapes. These are passed by name in JSON config and resolved at solve time.

## Git workflow

All development uses feature branches and pull requests. There is a single long-lived branch:

- **`master`** — stable, released code. All PRs target this branch.
- **feature branches** — cut from `master`, merged back via PR.

### Typical change

```bash
git checkout master
git pull
git checkout -b my-feature-branch
# ... make changes, commit ...
gh pr create --base master --title "..." --body "..."
```

### Release process

Before bumping, check the release against this list:

- Move the `[Unreleased]` entries in `CHANGELOG.md` and the "Unreleased" section of `docs/release-notes.qmd` under the new version number and date.
- Add the version to the docs of anything new in this release that was documented without one. For `clerq.autodiff` (first released in this version), add `.. versionadded:: <new version>` to the `AutodiffPropagation` class docstring and to `PropagationSolver.to_autodiff`, and "Requires clerq <new version> or later" to the top of `docs/examples/autodiff-pulse-design.ipynb`. After a docstring change, run `cd docs && uv run --extra docs quartodoc build` and commit the regenerated `docs/api`.

```bash
# On master: bump version, commit, tag, then push
uv run bump-my-version bump patch   # or minor / major
git push
git push --tags
```

The `git push --tags` triggers the CI publish workflow to PyPI.

## Environment

The primary workflow uses **uv** (see `pyproject.toml`). Run `uv sync --extra dev` to create the virtual environment and install all dependencies. Core runtime dependencies are QuTiP, NumPy, and SciPy. The package installs CLI entry points: `mbsolve` and `obsolve`.

An `environment.yml` is also kept as a conda fallback.
