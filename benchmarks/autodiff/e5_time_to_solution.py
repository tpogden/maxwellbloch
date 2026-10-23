"""E5: time to solution on a design problem: optimal storage in a Doppler-broadened Lambda medium.

The problem (see research/autodiff-notes.md): a probe is stored while the coupling field is off
(t = 3 to 5) and the energy retrieved after t = 5 is maximised. The medium is uniform with 11
Doppler classes (effective optical depth about 25). Two ways to design the fields:

  free-form   44 coefficients (20 for the coupling control, 24 for the probe shape), gradient
              descent (Adam) with exact gradients from reverse-mode autodiff
  baseline    5 parameters (a super-Gaussian probe and two coupling amplitudes), Nelder-Mead

and what the free-form design would cost with finite-difference gradients instead
(measured: the time of a 44-parameter forward-difference gradient).

Reports per-step times, the retrieved fraction against step, the step at which the free-form
design passes the baseline's best value, and the wall times. The longer runs behind the
numbers in the research notes (100 Adam steps) are not repeated in full: --steps sets the length
of the free-form run, --nm-evals that of the baseline search.
"""

import argparse
import json
import time

import numpy as np
from common import RESULTS, emit, summarise, time_once

T_OFF, T_ON, T_MIN = 3.0, 5.0, -3.0
CAP = 5.0
KC, KP = 20, 24


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def build():
    from clerq import PropagationSolver

    cfg = {
        "system": {
            "fields": [
                {"coupled_levels": [[0, 1]], "detuning": 0.0, "detuning_positive": True, "label": "probe",
                 "rabi_freq": 1.0e-3, "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": 1.0},
                 "rabi_freq_t_func": "gaussian"},
                {"coupled_levels": [[1, 2]], "detuning": 0.0, "detuning_positive": False, "label": "coupling",
                 "rabi_freq": CAP, "rabi_freq_t_args": {"ampl": 1.0, "fwhm": 0.2, "off": T_OFF, "on": T_ON},
                 "rabi_freq_t_func": "ramp_offon"},
            ],
            "num_states": 3,
            "decays": [{"channels": [[0, 1], [2, 1]], "rate": 1.0}, {"channels": [[2, 2]], "rate": 0.01}],
        },
        "t_min": T_MIN, "t_max": T_ON + 8.0, "t_steps": 140,
        "z_min": -0.2, "z_max": 1.2, "z_steps": 50, "z_steps_inner": 20,
        "interaction_strengths": [20.0, 20.0],
        "velocity_classes": {"thermal_delta_min": -4.0, "thermal_delta_max": 4.0,
                             "thermal_delta_steps": 9, "thermal_width": 2.0},
    }
    return PropagationSolver.from_json_str(json.dumps(cfg))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--steps", type=int, default=40)
    p.add_argument("--nm-evals", type=int, default=300)
    p.add_argument("--baseline-value", type=float, default=None,
                   help="the baseline's best retrieved fraction (default: the Nelder-Mead result of this run)")
    p.add_argument("--tag", default="e5")
    p.add_argument("--quick", action="store_true")
    a = p.parse_args()

    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    import optax
    from scipy.optimize import minimize

    from clerq.autodiff import controls

    out = RESULTS / f"{a.tag}.jsonl"
    solver = build()
    ad = solver.to_autodiff()
    tl = np.asarray(ad.tlist)
    omega0_default = np.asarray(ad.propagate())[:, 0, :].real
    E0 = float(np.sum(omega0_default[0] ** 2))
    hold = sigmoid((tl - T_OFF) / 0.1) * sigmoid((T_ON - tl) / 0.1)
    m_coupling, m_probe = 1 - hold, sigmoid((T_OFF - 0.6 - tl) / 0.2)
    Bc, Bp = controls.gaussian_basis(tl, KC), controls.gaussian_basis(tl, KP, hi=T_OFF)
    after = jnp.asarray(tl >= T_ON)

    def retrieved(om):
        O = ad.propagate({"omega0": om})
        return jnp.sum(after * jnp.abs(O[0, -1]) ** 2) / jnp.sum(jnp.abs(O[0, 0]) ** 2)

    def fields(coef):
        c = CAP * jax.nn.sigmoid(jnp.asarray(Bc) @ coef[:KC]) * jnp.asarray(m_coupling)
        pr = (jnp.asarray(Bp) @ coef[KC:]) * jnp.asarray(m_probe)
        return jnp.stack([pr * jnp.sqrt(E0 / jnp.sum(pr**2)), c])

    value = jax.jit(lambda coef: -retrieved(fields(coef)))
    vg = jax.jit(jax.value_and_grad(lambda coef: -retrieved(fields(coef))))
    rng = np.random.default_rng(1)
    coef0 = jnp.asarray(np.concatenate([rng.normal(0, 1.0, KC),
                                        np.linalg.lstsq(Bp, np.exp(-4 * np.log(2) * (tl / 1.2) ** 2), rcond=None)[0]
                                        * (1 + 0.3 * rng.normal(size=KP))]))
    n_par = KC + KP

    # --- cost of one 44-parameter gradient: reverse AD vs forward differences ---------------
    time_once(lambda: vg(coef0))
    t_ad = [time_once(lambda: vg(coef0)) for _ in range(3 if not a.quick else 1)]

    def fd_grad():
        f0 = value(coef0)
        return jnp.stack([(value(coef0 + 1e-6 * jnp.eye(n_par)[i]) - f0) / 1e-6 for i in range(n_par)])

    time_once(fd_grad)
    t_fd = [time_once(fd_grad) for _ in range(2 if not a.quick else 1)]
    t_fwd = [time_once(lambda: value(coef0)) for _ in range(3)]
    emit({"exp": a.tag, "what": "gradient_cost", "n_par": n_par, "forward": summarise(t_fwd),
          "reverse_ad": summarise(t_ad), "fd_forward_differences": summarise(t_fd)}, out)
    print("44-parameter gradient: AD", round(np.median(t_ad), 2), "s, FD", round(np.median(t_fd), 2), "s", flush=True)

    # --- baseline: 5-parameter family, Nelder-Mead ------------------------------------------
    def fam(x):
        w, c, aw, ar, kap = abs(x[0]) + 0.2, x[1], np.clip(x[2], 0, CAP), np.clip(x[3], 0, CAP), np.clip(x[4], 1, 8)
        pr = np.exp(-np.log(2) * (2 * np.abs(tl - c) / w) ** kap) * m_probe
        pr = pr * np.sqrt(E0 / np.sum(pr**2))
        ramp = sigmoid((tl - T_ON) / 0.1)
        return jnp.asarray(np.stack([pr, (aw * (1 - ramp) + ar * ramp) * m_coupling]))

    retr = jax.jit(retrieved)
    nm_history = []
    t_nm0 = time.perf_counter()

    def nm_obj(x):
        v = -float(retr(fam(x)))
        nm_history.append((time.perf_counter() - t_nm0, -v))
        return v

    # a generic start, not the optimum found in earlier runs
    res = minimize(nm_obj, [1.0, 2.4, 4.0, 3.0, 1.5], method="Nelder-Mead",
                   options=dict(maxfev=10 if a.quick else a.nm_evals, xatol=1e-3, fatol=1e-6))
    best_nm = -res.fun
    t_nm = time.perf_counter() - t_nm0
    emit({"exp": a.tag, "what": "baseline_nelder_mead", "n_par": 5, "evals": int(res.nfev), "wall_s": t_nm,
          "retrieved": best_nm, "history": nm_history[:: max(1, len(nm_history) // 50)]}, out)
    print("baseline: retrieved", round(best_nm, 4), "after", res.nfev, "evals,", round(t_nm), "s", flush=True)
    target = a.baseline_value if a.baseline_value is not None else best_nm

    # --- free-form design: Adam with exact gradients -----------------------------------------
    steps = 3 if a.quick else a.steps
    opt = optax.adam(optax.cosine_decay_schedule(0.1, max(steps, 100), alpha=0.05))
    state, coef = opt.init(coef0), coef0
    hist, t0 = [], time.perf_counter()
    for k in range(steps):
        v, g = vg(coef)
        hist.append((time.perf_counter() - t0, -float(v)))
        upd, state = opt.update(g, state)
        coef = optax.apply_updates(coef, upd)
    final = -float(value(coef))
    t_adam = time.perf_counter() - t0
    first_pass = next((k for k, (_, r) in enumerate(hist) if r >= target), None)
    emit({"exp": a.tag, "what": "free_form_adam", "n_par": n_par, "steps": steps, "wall_s": t_adam,
          "retrieved_final": final, "history": hist, "target": target, "step_passing_target": first_pass,
          "wall_s_to_target": None if first_pass is None else hist[first_pass][0],
          "fd_gradient_wall_s_to_target_estimate": None if first_pass is None else (first_pass + 1) * float(np.median(t_fd))}, out)
    print("Adam:", steps, "steps", round(t_adam), "s, final", round(final, 4), "passes the baseline at step", first_pass, flush=True)


if __name__ == "__main__":
    main()
