# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

MaxwellBloch is a Python package for numerically solving the coupled Maxwell-Bloch equations, which describe nonlinear propagation of near-resonant light through thermal atomic vapors. It models two, three, and many-level quantum systems using density matrix formalism via the QuTiP library.

## Commands

```bash
# Run all tests in parallel
pytest -n auto

# Run a single test file
pytest maxwellbloch/tests/test_mb_solve.py

# Run a single test
pytest maxwellbloch/tests/test_mb_solve.py::TestClassName::test_method_name

# Run with coverage
pytest --cov -n auto

# Build docs
sphinx-build docs docs/_build -b html

# Build distribution
python setup.py sdist --formats=gztar bdist_wheel
```

Make targets: `test`, `test_cov`, `docs_html`, `dist`, `clean_qu` (removes cached `.qu` files).

## Architecture

### Class hierarchy

```
OBBase (ob_base.py)          — QuTiP interface, operators, density matrix math
  └── OBSolve (ob_solve.py)  — time-domain master equation solver (wraps qutip.mesolve)
        └── MBSolve (mb_solve.py) — spatial propagation through a medium

OBAtom (ob_atom.py)          — atomic system: energy levels, decay channels, coupled fields
Field (field.py)             — a laser field coupling two atomic levels (Rabi freq, detuning)
```

Hyperfine structure lives in `hyperfine.py` and `angmom.py` (3j/6j symbols, Clebsch-Gordan coefficients).

### Solver flow

- **OBSolve** calls `qutip.mesolve` to evolve the atomic density matrix in time for a single spatial point.
- **MBSolve** wraps OBSolve and propagates the field through space using finite differences (Adams or Euler method). At each z-step it: (1) solves the Bloch equations, (2) reads off coherences, (3) updates the Rabi frequency for the next step.
- Velocity classes are used for Doppler broadening: OBSolve is run independently per velocity class and results are summed.

### Configuration

Problems are defined as JSON (see `fileio.py`). `MBSolve` and `OBSolve` accept a JSON string or a dict. Test fixtures live in `maxwellbloch/tests/json/`.

### Caching

Solved results are pickled by QuTiP into `.qu` files next to the script. Use `make clean_qu` to remove them.

### Time functions

`t_funcs.py` provides named callables (sech, Gaussian, square) for pulse shapes. These are passed by name in JSON config and resolved at solve time.

## Environment

Use the conda environment defined in `environment.yml`. Core runtime dependencies are QuTiP, NumPy, and SciPy. The package installs CLI entry points: `mbsolve` and `obsolve`.
