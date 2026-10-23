"""Shared pieces of the ``clerq.autodiff`` benchmarks: systems, timing, statistics, records.

The benchmark problems all use a free-form probe control, ``omega0[0] = B @ coef`` with
``B`` a Gaussian basis from :mod:`clerq.autodiff.controls`, so the number of parameters
``N`` is the number of coefficients and can be varied freely. Other fields keep their
analytic input.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

RESULTS = Path(__file__).parent / "results"


def _run(cmd: list[str]) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:  # pragma: no cover - metadata only
        return ""


def meta() -> dict:
    """Hardware, software and machine state to store with every record."""
    import diffrax
    import jax

    import clerq

    src = Path(clerq.__file__).resolve().parent
    return {
        "date": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "git_sha": _run(["git", "-C", str(src), "rev-parse", "--short", "HEAD"]),
        "clerq_src": str(src),
        "clerq": getattr(clerq, "__version__", ""),
        "python": sys.version.split()[0],
        "jax": jax.__version__,
        "diffrax": diffrax.__version__,
        "devices": [str(d) for d in jax.devices()],
        "machine": platform.platform(),
        "cpu": _run(["sysctl", "-n", "machdep.cpu.brand_string"]),
        "cpu_count": os.cpu_count(),
        "load_avg": [round(x, 2) for x in os.getloadavg()],
        "power": (_run(["pmset", "-g", "batt"]).splitlines() or [""])[0],
    }


def emit(record: dict, path: str | Path) -> None:
    """Append one JSON record (with metadata) to ``path``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps({**record, "meta": meta()}) + "\n")


def block(x):
    import jax

    return jax.block_until_ready(x)


def time_once(fn) -> float:
    t = time.perf_counter()
    block(fn())
    return time.perf_counter() - t


def summarise(times: list[float], seed: int = 0) -> dict:
    """Median and a 95% bootstrap interval of the median (min/max too)."""
    a = np.asarray(times, dtype=float)
    out = {"n": len(a), "median": float(np.median(a)), "min": float(a.min()), "max": float(a.max())}
    if len(a) >= 5:
        rng = np.random.default_rng(seed)
        boots = np.median(rng.choice(a, size=(2000, len(a)), replace=True), axis=1)
        out["ci_lo"], out["ci_hi"] = (float(x) for x in np.percentile(boots, [2.5, 97.5]))
    return out


# ---- systems ------------------------------------------------------------------------


def _solver(cfg: dict):
    from clerq import PropagationSolver

    return PropagationSolver.from_json_str(json.dumps(cfg))


def two_level(n_t: int = 1000, n_z: int = 200):
    """A resonant sech 2-pi pulse through a two-level medium."""
    return _solver(
        {
            "system": {
                "num_states": 2,
                "fields": [
                    {
                        "coupled_levels": [[0, 1]],
                        "rabi_freq_t_func": "sech",
                        "rabi_freq_t_args": {"ampl": 0.8384014365421667, "centre": 0.0, "width": 0.3796628587572578},
                    }
                ],
            },
            "t_min": -2.0, "t_max": 10.0, "t_steps": n_t,
            "z_min": 0.0, "z_max": 1.0, "z_steps": n_z, "z_steps_inner": 2,
            "interaction_strengths": [10.0],
        }
    )


def lambda_eit(n_t: int = 200, n_z: int = 60):
    """A weak Gaussian probe and a ramped coupling field in a lambda medium."""
    return _solver(
        {
            "system": {
                "num_states": 3,
                "decays": [{"channels": [[0, 1], [1, 2]], "rate": 1.0}],
                "fields": [
                    {
                        "coupled_levels": [[0, 1]],
                        "rabi_freq_t_func": "gaussian",
                        "rabi_freq_t_args": {"ampl": 0.1, "centre": 0.0, "fwhm": 0.1},
                    },
                    {
                        "coupled_levels": [[1, 2]],
                        "rabi_freq_t_func": "ramp_offon",
                        "rabi_freq_t_args": {"ampl": 10.0, "fwhm": 0.05, "off": 0.3, "on": 0.7},
                    },
                ],
            },
            "t_min": 0.0, "t_max": 1.0, "t_steps": n_t,
            "z_min": -0.2, "z_max": 1.2, "z_steps": n_z, "z_steps_inner": 2,
            "interaction_strengths": [20.0, 20.0],
        }
    )


def rb87(n_t: int = 200, n_z: int = 50, n_classes_steps: int = 9):
    """Rb87 5S1/2 F=1,2 and 5P1/2 F=1 (11 sublevels) with Doppler classes (11 by default)."""
    from clerq import hyperfine

    atom = hyperfine.Atom1e(element="Rb", isotope="87")
    atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))
    atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=2))
    atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))
    cfg = {
        "system": {
            "num_states": atom.get_num_mF_levels(),
            "energies": atom.get_energies(),
            "decays": [
                {
                    "rate": 1.0,
                    "channels": atom.get_coupled_levels(F_level_idxs_a=(0, 1), F_level_idxs_b=(2,)),
                    "factors": atom.get_decay_factors(F_level_idxs_a=(0, 1), F_level_idxs_b=(2,)).tolist(),
                }
            ],
            "fields": [
                {
                    "coupled_levels": atom.get_coupled_levels(F_level_idxs_a=(0,), F_level_idxs_b=(2,)),
                    "factors": atom.get_clebsch_hf_factors(F_level_idxs_a=(0,), F_level_idxs_b=(2,), q=1).tolist(),
                    "detuning": 0.0,
                    "rabi_freq": 1e-3,
                    "rabi_freq_t_func": "gaussian",
                    "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": 1.0},
                }
            ],
            "initial_state": [0.5 / 3.0] * 3 + [0.5 / 5.0] * 5 + [0.0] * 3,
        },
        "t_min": -2.0, "t_max": 6.0, "t_steps": n_t,
        "z_min": 0.0, "z_max": 1.0, "z_steps": n_z, "z_steps_inner": 2,
        "interaction_strengths": [14.0],
    }
    if n_classes_steps:
        cfg["velocity_classes"] = {
            "thermal_delta_min": -3.0, "thermal_delta_max": 3.0,
            "thermal_delta_steps": n_classes_steps, "thermal_width": 1.5,
        }
    return _solver(cfg)


SYSTEMS = {"two_level": two_level, "lambda_eit": lambda_eit, "rb87": rb87}


def gaussian_basis(t, n: int, lo=None, hi=None):
    """Gaussian bumps on a regular grid of centres (a copy of clerq.autodiff.controls.gaussian_basis,
    so that the harness also runs against source trees that predate it)."""
    t = np.asarray(t, dtype=float)
    lo = t[0] if lo is None else lo
    hi = t[-1] if hi is None else hi
    centres = np.linspace(lo, hi, n)
    sigma = centres[1] - centres[0]
    return np.exp(-0.5 * ((t[:, None] - centres[None, :]) / sigma) ** 2)


class Problem:
    """A benchmark problem: a solver, its autodiff object and a coefficient-parametrised loss."""

    def __init__(self, solver, n_coef: int, **ad_options):
        import jax
        import jax.numpy as jnp

        self.solver = solver
        self.ad = solver.to_autodiff(**ad_options)
        self.n_coef = n_coef
        t = np.asarray(self.ad.tlist)
        first = np.asarray(self.ad.propagate())[:, 0, :].real  # analytic inputs, cyclic
        self.om_fixed = jnp.asarray(first)
        self.B = gaussian_basis(t, n_coef)
        self.coef0 = jnp.asarray(np.linalg.lstsq(self.B, first[0], rcond=None)[0])
        B = jnp.asarray(self.B)
        om_fixed = self.om_fixed

        def field(coef):
            return om_fixed.at[0].set(B @ coef)

        self._field = field

        def loss_with(propagate):
            def loss(coef):
                O = propagate({"omega0": field(coef)})
                return jnp.sum(jnp.abs(O[0, -1]) ** 2)

            return loss

        self.loss = loss_with(self.ad.propagate)
        self.value_and_grad = jax.jit(jax.value_and_grad(self.loss))
        self.value = jax.jit(self.loss)

    def forward_mode_grad(self):
        import diffrax as dfx
        import jax
        import jax.numpy as jnp

        ad, field = self.ad, self._field

        def loss(coef):
            O = ad._propagate({"omega0": field(coef)}, adjoint=dfx.ForwardMode())
            return jnp.sum(jnp.abs(O[0, -1]) ** 2)

        return jax.jit(lambda c: (loss(c), jax.jacfwd(loss)(c)))

    def fd_jax(self):
        """Forward-difference gradient through the compiled loss (N + 1 evaluations)."""
        import jax.numpy as jnp

        value, n, h = self.value, self.n_coef, 1e-6

        def grad(coef):
            f0 = value(coef)
            return jnp.stack([(value(coef + h * jnp.eye(n)[i]) - f0) / h for i in range(n)])

        return grad


def make_solver(name: str, quick: bool = False):
    if quick:
        small = {"two_level": (100, 20), "lambda_eit": (60, 10), "rb87": (60, 10)}[name]
        return SYSTEMS[name](*small)
    return SYSTEMS[name]()
