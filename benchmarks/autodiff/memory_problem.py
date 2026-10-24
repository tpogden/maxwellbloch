"""The quantum-memory design problem used by E5 and E7.

A probe pulse is stored in a uniform lambda medium while the coupling field is off (t = 3 to 5), and
the energy retrieved after t = 5 is maximised. The coupling control (20 coefficients) and the probe
shape (24 coefficients) are expanded in Gaussian bumps: 44 parameters. See research/autodiff-notes.md.
"""

import json

import numpy as np

T_OFF, T_ON, T_MIN = 3.0, 5.0, -3.0
CAP = 5.0
KC, KP = 20, 24


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def solver(g: float = 12.5, doppler: bool = False):
    """The medium. ``g`` = 12.5 gives an optical depth of 25 without Doppler classes."""
    from clerq import PropagationSolver

    cfg = {
        "system": {
            "fields": [
                {
                    "coupled_levels": [[0, 1]],
                    "detuning": 0.0,
                    "detuning_positive": True,
                    "label": "probe",
                    "rabi_freq": 1.0e-3,
                    "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": 1.0},
                    "rabi_freq_t_func": "gaussian",
                },
                {
                    "coupled_levels": [[1, 2]],
                    "detuning": 0.0,
                    "detuning_positive": False,
                    "label": "coupling",
                    "rabi_freq": CAP,
                    "rabi_freq_t_args": {
                        "ampl": 1.0,
                        "fwhm": 0.2,
                        "off": T_OFF,
                        "on": T_ON,
                    },
                    "rabi_freq_t_func": "ramp_offon",
                },
            ],
            "num_states": 3,
            "decays": [
                {"channels": [[0, 1], [2, 1]], "rate": 1.0},
                {"channels": [[2, 2]], "rate": 0.01},
            ],
        },
        "t_min": T_MIN,
        "t_max": T_ON + 8.0,
        "t_steps": 140,
        "z_min": -0.2,
        "z_max": 1.2,
        "z_steps": 50,
        "z_steps_inner": 20,
        "interaction_strengths": [g, g],
    }
    if doppler:
        cfg["velocity_classes"] = {
            "thermal_delta_min": -4.0,
            "thermal_delta_max": 4.0,
            "thermal_delta_steps": 9,
            "thermal_width": 2.0,
        }
    return PropagationSolver.from_json_str(json.dumps(cfg))


class MemoryProblem:
    """Objective ``-retrieved(coef)`` of the 44 coefficients, with a jitted value and gradient."""

    def __init__(self, g: float = 12.5, doppler: bool = False):
        import jax
        import jax.numpy as jnp

        self.ad = solver(g, doppler).to_autodiff()
        tl = np.asarray(self.ad.tlist)
        self.tl = tl
        omega0_default = np.asarray(self.ad.propagate())[:, 0, :].real
        E0 = float(np.sum(omega0_default[0] ** 2))
        hold = sigmoid((tl - T_OFF) / 0.1) * sigmoid((T_ON - tl) / 0.1)
        m_coupling, m_probe = (
            jnp.asarray(1 - hold),
            jnp.asarray(sigmoid((T_OFF - 0.6 - tl) / 0.2)),
        )
        from common import gaussian_basis

        self.Bc = gaussian_basis(tl, KC)
        self.Bp = gaussian_basis(tl, KP, hi=T_OFF)
        Bc, Bp = jnp.asarray(self.Bc), jnp.asarray(self.Bp)
        after = jnp.asarray(tl >= T_ON)
        ad = self.ad

        def retrieved(om):
            field_zt = ad.propagate({"omega0": om})
            return jnp.sum(after * jnp.abs(field_zt[0, -1]) ** 2) / jnp.sum(
                jnp.abs(field_zt[0, 0]) ** 2
            )

        def objective(coef):
            c = CAP * jax.nn.sigmoid(Bc @ coef[:KC]) * m_coupling
            pr = (Bp @ coef[KC:]) * m_probe
            return -retrieved(jnp.stack([pr * jnp.sqrt(E0 / jnp.sum(pr**2)), c]))

        self.value = jax.jit(objective)
        self.value_and_grad = jax.jit(jax.value_and_grad(objective))
        self.n_par = KC + KP

    def random_start(self, seed: int) -> np.ndarray:
        """Coupling logits N(0, 1) and a noisy Gaussian probe, as in the E5 free-form run."""
        rng = np.random.default_rng(seed)
        probe = np.linalg.lstsq(
            self.Bp, np.exp(-4 * np.log(2) * (self.tl / 1.2) ** 2), rcond=None
        )[0]
        return np.concatenate(
            [rng.normal(0, 1.0, KC), probe * (1 + 0.3 * rng.normal(size=KP))]
        )
