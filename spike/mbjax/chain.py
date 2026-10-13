"""A parameterised propagation problem for benchmarking: builds ``loss(theta)``.

theta = [g_scale, gamma_scale, detuning, c_1..c_4]: interaction strength and decay-rate
scales, a fixed-detuning shift, and four Gaussian-basis coefficients of the input pulse.
Loss = integral |Omega_out/2pi - target|^2 dt (target = delayed, attenuated pulse).
"""

from __future__ import annotations

import diffrax as dfx
import jax
import jax.numpy as jnp
import numpy as np

from mbjax.lindblad import gaussian
from mbjax.propagate import build_propagation_model, propagate, velocity_classes
from mbjax.systems import make_system

jax.config.update("jax_enable_x64", True)

# (pulse ampl in cyclic units, interaction strength g): chosen so the medium visibly acts
SYSTEM_PARAMS = {"two_level": (1.0, 0.5), "j12_12": (1.5, 0.4), "j1_2_res": (2.5, 0.4),
                 "j72_72": (3.0, 1.5)}
T_MIN, T_MAX = -3.0, 9.0
BASIS_CENTRES = np.array([0.0, 0.6, 1.2, 1.8])


def make_chain(system: str, n_z: int, n_t: int, n_v: int = 1, *, strategy: str = "s3",
               rtol: float = 1e-6, atol: float = 1e-8, max_steps: int = 16384,
               checkpoints: int | None = None, chunk: int | None = None,
               forward_mode: bool = False, n_steps: int | None = None, solver=None):
    ampl, g = SYSTEM_PARAMS[system]
    spec = make_system(system)
    spec["fields"][0].update(rabi_freq=1.0, rabi_freq_t_func="gaussian",
                             rabi_freq_t_args={"ampl": ampl, "centre": 0.0, "fwhm": 1.0})
    model = build_propagation_model(spec)
    ts = jnp.linspace(T_MIN, T_MAX, n_t + 1)
    h = 1.0 / n_z
    N_z = jnp.ones(n_z)
    vc = None if n_v == 1 else {"thermal_delta_min": -1.0, "thermal_delta_max": 1.0,
                                "thermal_delta_steps": n_v - 1, "thermal_width": 0.7,
                                # clerq adds an inner grid (default one point at 0); put it
                                # on an outer point so exactly n_v classes remain
                                "thermal_delta_inner_min": -1.0, "thermal_delta_inner_max": -1.0}
    deltas, w = velocity_classes(vc)
    basis = jnp.stack([ampl * gaussian(1.0, c, 1.0)(ts) for c in BASIS_CENTRES])
    target = 0.5 * ampl * gaussian(1.0, 1.0, 1.2)(ts)
    if forward_mode:
        adjoint = dfx.ForwardMode()
    elif checkpoints is not None:
        adjoint = dfx.RecursiveCheckpointAdjoint(checkpoints=checkpoints)
    else:
        adjoint = None

    def loss(theta):
        coef = theta[3:]
        om0 = 2 * jnp.pi * (coef @ basis)
        m = dict(model, detuning=theta[2])
        om = propagate(m, om0, ts, N_z, g * theta[0], h, deltas, w,
                       strategy="s2" if strategy in ("s1", "fwd") else strategy,
                       gamma_scale=theta[1], rtol=rtol, atol=atol, n_steps=n_steps,
                       adjoint=adjoint, chunk=chunk, max_steps=max_steps, solver=solver)
        out = om[-1].real / (2 * jnp.pi)
        return jnp.sum((out - target) ** 2) * ((T_MAX - T_MIN) / n_t)

    theta0 = jnp.array([1.0, 1.0, 0.0, 1.0, 0.1, 0.0, 0.0])
    return loss, theta0
