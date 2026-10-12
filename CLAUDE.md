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

### Release notes

`docs/release-notes.qmd` is the user-readable version of `CHANGELOG.md`. When a PR adds a CHANGELOG entry, update the
matching section of the release notes (the "Unreleased" section until a release is cut, then retitle it with the version
and date).

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

### Worktrees

Do all code changes in a dedicated git worktree cut from `next`, not in the main checkout. At the start of a task, use EnterWorktree with a branch name describing the change. Run `uv sync --extra dev` in the new worktree before running tests, since `.venv` and `.qu` caches are not shared.

## Environment

The primary workflow uses **uv** (see `pyproject.toml`). Run `uv sync --extra dev` to create the virtual environment and install all dependencies. Core runtime dependencies are QuTiP, NumPy, and SciPy. The package installs CLI entry points: `mbsolve` and `obsolve`.

An `environment.yml` is also kept as a conda fallback.
