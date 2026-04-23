# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

MaxwellBloch is a Python package for numerically solving the coupled Maxwell-Bloch equations, which describe nonlinear propagation of near-resonant light through thermal atomic vapors. It models two, three, and many-level quantum systems using density matrix formalism via the QuTiP library.

## Commands

```bash
# Create/sync environment (first time or after dependency changes)
uv sync --extra dev

# Run all tests in parallel
uv run pytest -n auto

# Run a single test file
uv run pytest tests/test_mb_solve.py

# Run a single test
uv run pytest tests/test_mb_solve.py::TestClassName::test_method_name

# Run with coverage
uv run pytest --cov -n auto

# Run performance benchmarks
uv run pytest tests/bench_mb_solve.py --benchmark-only

# Lint
uv run ruff check .

# Format
uv run ruff format .

# Build docs
uv run sphinx-build docs docs/_build -b html

# Build distribution
uv build

# Release — bump version, commit, tag, then push to trigger CI publish
uv run bump-my-version bump patch   # bug fixes
uv run bump-my-version bump minor   # new features, backwards-compatible
uv run bump-my-version bump major   # breaking changes
git push && git push --tags
```

Make targets: `test`, `test_cov`, `bench`, `lint`, `format`, `format_check`, `docs_html`, `dist`, `bump_patch`, `bump_minor`, `bump_major`, `clean_qu` (removes cached `.qu` files).

## Architecture

### Class hierarchy

```
OBBase (ob_base.py)          — QuTiP interface, operators, density matrix math
  └── OBAtom (ob_atom.py)    — atomic system: energy levels, decay channels, coupled fields

OBSolve (ob_solve.py)        — time-domain master equation solver (wraps qutip.mesolve)
  └── MBSolve (mb_solve.py)  — spatial propagation through a medium

Field (field.py)             — a laser field coupling two atomic levels (Rabi freq, detuning)
```

OBSolve holds an `OBAtom` instance as `self.atom` and delegates Hamiltonian construction to it.

Hyperfine structure lives in `hyperfine.py` and `angmom.py` (3j/6j symbols, Clebsch-Gordan coefficients).

### Solver flow

- **OBSolve** calls `qutip.mesolve` to evolve the atomic density matrix in time for a single spatial point.
- **MBSolve** wraps OBSolve and propagates the field through space using finite differences (Adams or Euler method). At each z-step it: (1) solves the Bloch equations, (2) reads off coherences, (3) updates the Rabi frequency for the next step.
- Velocity classes are used for Doppler broadening: OBSolve is run independently per velocity class and results are summed.

### Configuration

Problems are defined as JSON (see `fileio.py`). `MBSolve` and `OBSolve` accept a JSON string or a dict. Test fixtures live in `tests/json/`.

### Caching

Solved results are pickled by QuTiP into `.qu` files next to the script. Use `make clean_qu` to remove them.

### Time functions

`t_funcs.py` provides named callables (sech, Gaussian, square) for pulse shapes. These are passed by name in JSON config and resolved at solve time.

## Git workflow

All development uses feature branches and pull requests. The branch structure is:

- **`master`** — stable, released code. Only receives merges from `next` as part of a release.
- **`next`** — the release branch. All PRs target this branch.
- **feature branches** — cut from `next`, merged back via PR.

### Typical change

```bash
git checkout next
git pull
git checkout -b my-feature-branch
# ... make changes, commit ...
gh pr create --base next --title "..." --body "..."
```

### Release process

```bash
# On next: bump version, commit, then merge to master and tag
uv run bump-my-version bump patch   # or minor / major
git checkout master
git merge --no-ff next
git push
git push --tags
git checkout next
```

The `git push --tags` triggers the CI publish workflow to PyPI.

## Environment

The primary workflow uses **uv** (see `pyproject.toml`). Run `uv sync --extra dev` to create the virtual environment and install all dependencies. Core runtime dependencies are QuTiP, NumPy, and SciPy. The package installs CLI entry points: `mbsolve` and `obsolve`.

An `environment.yml` is also kept as a conda fallback.
