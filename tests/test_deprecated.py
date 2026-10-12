"""Tests that the names from before the API rename still work, with a
DeprecationWarning."""

import json
import os

import numpy as np
import pytest

import clerq
from clerq import LindbladSolver, PropagationSolver, System, system_base

JSON_DIR = os.path.abspath(os.path.join(__file__, "../", "json"))

SYSTEM = {
    "num_states": 2,
    "decays": [{"channels": [[0, 1]], "rate": 1.0}],
    "fields": [{"coupled_levels": [[0, 1]], "rabi_freq": 1.0e-3}],
}


@pytest.mark.parametrize(
    ("module", "name", "new"),
    [
        ("mb_solve", "MBSolve", PropagationSolver),
        ("ob_solve", "OBSolve", LindbladSolver),
        ("ob_atom", "OBAtom", System),
        ("ob_base", "OBBase", system_base.SystemBase),
    ],
)
def test_old_class_names(module, name, new):
    """The old modules still import and the old class names are aliases."""
    old_module = __import__(f"clerq.{module}", fromlist=[name])
    with pytest.warns(DeprecationWarning, match=name):
        assert getattr(old_module, name) is new
    with pytest.raises(AttributeError):
        old_module.not_a_name  # noqa: B018


def test_top_level_exports():
    assert clerq.__all__ == ["Field", "LindbladSolver", "PropagationSolver", "System"]
    for name in clerq.__all__:
        assert getattr(clerq, name).__name__ == name


@pytest.mark.parametrize("cls", [LindbladSolver, PropagationSolver])
class TestAtomIsSystem:
    def test_atom_arg(self, cls):
        with pytest.warns(DeprecationWarning, match="'atom' argument"):
            solver = cls(atom=SYSTEM)
        assert solver.system.num_states == 2

    def test_atom_and_system_args(self, cls):
        with pytest.raises(TypeError):
            cls(system=SYSTEM, atom=SYSTEM)

    def test_atom_attribute(self, cls):
        solver = cls(system=SYSTEM)
        with pytest.warns(DeprecationWarning, match="'atom' attribute"):
            assert solver.atom is solver.system

    def test_atom_key_in_json_str(self, cls):
        with pytest.warns(DeprecationWarning, match="'atom' key"):
            solver = cls.from_json_str(json.dumps({"atom": SYSTEM}))
        assert solver.system.num_states == 2
        assert "system" in solver.get_json_dict()
        assert "atom" not in solver.get_json_dict()

    def test_atom_key_in_json_file(self, cls, tmp_path):
        path = tmp_path / "problem.json"
        path.write_text(json.dumps({"atom": SYSTEM}))
        with pytest.warns(DeprecationWarning, match="'atom' key"):
            solver = cls.from_json(str(path))
        assert solver.system.num_states == 2

    def test_atom_and_system_keys(self, cls):
        with pytest.raises(ValueError):
            cls.from_json_str(json.dumps({"atom": SYSTEM, "system": SYSTEM}))


def test_obsolve():
    solver = LindbladSolver(system=SYSTEM, t_steps=10)
    states_t = solver.solve()
    with pytest.warns(DeprecationWarning, match="obsolve"):
        assert np.array_equal(solver.obsolve(), states_t)


def test_obsolve_on_propagation_solver():
    """obsolve() on a PropagationSolver solved a single point; it still does."""
    solver = PropagationSolver(system=SYSTEM, t_steps=10, interaction_strengths=[1.0])
    with pytest.warns(DeprecationWarning, match="obsolve"):
        states_t = solver.obsolve()
    assert states_t.shape == (11, 2, 2)


@pytest.mark.parametrize(
    ("old", "step"),
    [("mbsolve", "ab"), ("mbsolve_ab", "ab"), ("mbsolve_euler", "euler")],
)
def test_mbsolve(old, step):
    kwargs = {"system": SYSTEM, "t_steps": 10, "z_steps": 4}
    solver = PropagationSolver(interaction_strengths=[1.0], **kwargs)
    Omegas_zt, states_zt = solver.solve(step=step, progress=False)
    old_solver = PropagationSolver(interaction_strengths=[1.0], **kwargs)
    with pytest.warns(DeprecationWarning, match=old):
        old_Omegas_zt, old_states_zt = getattr(old_solver, old)(progress=False)
    assert np.array_equal(old_Omegas_zt, Omegas_zt)
    assert np.array_equal(old_states_zt, states_zt)


class TestSolverArg:
    """Functions that take a solver named the argument mb_solve, mbs or
    mbs_list."""

    @pytest.fixture(scope="class")
    def solver(self):
        solver = PropagationSolver(
            system=SYSTEM, t_steps=10, z_steps=4, interaction_strengths=[1.0]
        )
        solver.solve(progress=False)
        return solver

    def test_mb_solve(self, solver):
        from clerq import fixed, spectral

        for func, kwargs in [
            (spectral.freq_list, {}),
            (spectral.rabi_freq, {"field_idx": 0}),
            (spectral.absorption, {"field_idx": 0}),
            (spectral.dispersion, {"field_idx": 0}),
            (fixed.t_list, {"speed_of_light": 10.0}),
            (fixed.rabi_freq, {"field_idx": 0, "speed_of_light": 10.0}),
        ]:
            expected = func(solver, **kwargs)
            assert np.array_equal(func(solver=solver, **kwargs), expected)
            with pytest.warns(DeprecationWarning, match="'mb_solve' argument"):
                assert np.array_equal(func(mb_solve=solver, **kwargs), expected)

    def test_mb_solve_fileio(self, solver, tmp_path):
        from clerq import fileio

        for func in [fileio.save_csv_rabi_freq, fileio.save_csv_rabi_freq_abs]:
            path = tmp_path / f"{func.__name__}.csv"
            with pytest.warns(DeprecationWarning, match="'mb_solve' argument"):
                func(mb_solve=solver, filename=str(path))
            assert path.exists()

    def test_mbs(self, solver):
        from clerq import plot
        from clerq.plot import utils

        for func in [
            plot.field_spacetime,
            plot.field_envelope,
            plot.pulse_area,
            plot.population,
            plot.spectrum,
            utils.omega_abs,
            utils.pulse_area_z,
            utils.field_label,
        ]:
            func(solver=solver)
            with pytest.warns(DeprecationWarning, match="'mbs' argument"):
                func(mbs=solver)

    def test_mbs_list(self, solver):
        from clerq import plot

        plot.spectrum_overlay(solvers=[solver])
        with pytest.warns(DeprecationWarning, match="'mbs_list' argument"):
            plot.spectrum_overlay(mbs_list=[solver])

    def test_old_and_new(self, solver):
        from clerq import spectral

        with pytest.raises(TypeError):
            spectral.freq_list(solver=solver, mb_solve=solver)
