"""T0.4: reference problem from clerq (PropagationSolver) + convergence study.

Reference problem (two-level, Gaussian probe, square medium on z in [0, 1]):
    Gamma = 1.0, fwhm = 1.0, peak Rabi = 2*pi*ampl (cyclic ampl), t in [t_min, t_max]
Conventions pinned here (see PLAN.md T0.4):
    * Omegas_zt is the complex Rabi frequency in angular units; only Re(Omega) drives H.
    * Euler z-step: Omega_{j+1} = Omega_j + h * 2i * N(z_{j+1}) * 2*pi*g * sum_coh(z_j).
    * Adams-Bashforth bootstraps with one Euler step.
    * states_zt[j+1] holds the state solved at z_j (off-by-one); only Omegas_zt is used.
"""

import json
import sys
import time
import warnings

import numpy as np
from clerq import PropagationSolver

sys.path.insert(0, ".")
from mbjax.systems import make_system  # noqa: E402

warnings.filterwarnings("ignore")

PROBLEM = dict(
    system="two_level", gamma=1.0, ampl=1.0, fwhm=1.0, t_min=-3.0, t_max=9.0, g=0.5
)


def run(z_steps=100, inner=1, t_steps=1000, step="ab", rtol=1e-6, atol=1e-8,
        problem=None):
    p = {**PROBLEM, **(problem or {})}
    spec = make_system(p["system"], gamma=p["gamma"])
    spec["fields"][0].update(
        rabi_freq=1.0,
        rabi_freq_t_func="gaussian",
        rabi_freq_t_args={"ampl": p["ampl"], "centre": 0.0, "fwhm": p["fwhm"]},
    )
    s = PropagationSolver(
        system=spec, t_min=p["t_min"], t_max=p["t_max"], t_steps=t_steps,
        z_min=0.0, z_max=1.0, z_steps=z_steps, z_steps_inner=inner,
        interaction_strengths=[p["g"]],
        opts={"atol": atol, "rtol": rtol, "method": "adams", "nsteps": 100000},
    )
    t0 = time.perf_counter()
    s.solve(step=step, recalc=True, progress=False)
    return s.tlist, s.zlist, s.Omegas_zt[0].real.copy(), time.perf_counter() - t0


def rel_l2(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def sub(arr_t, n_full, n):
    """Indices selecting an n-step grid from an n_full-step grid (n divides n_full)."""
    return np.arange(0, n_full + 1, n_full // n)


def convergence():
    T_TRUTH, Z_TRUTH = 2000, 3200
    out = {"problem": PROBLEM, "truth": {"z_steps": Z_TRUTH, "t_steps": T_TRUTH,
                                         "step": "ab", "rtol": 1e-10, "atol": 1e-12}}
    t, z, om, dt = run(z_steps=Z_TRUTH, t_steps=T_TRUTH, step="ab", rtol=1e-10, atol=1e-12)
    truth = om[-1]
    out["truth"]["seconds"] = dt
    np.savez("refs/two_level_truth.npz", t=t, z=z, omega_zt=om, **{"problem": json.dumps(PROBLEM)})
    print("truth", dt, "s", flush=True)

    def err(tt, om_, nt):
        idx = sub(tt, T_TRUTH, nt)
        return rel_l2(om_[-1], truth[idx])

    A = []
    for step in ("euler", "ab"):
        for zs in (25, 50, 100, 200, 400, 800, 1600):
            tt, zz, o, dt = run(z_steps=zs, t_steps=1000, step=step, rtol=1e-8, atol=1e-10)
            A.append({"step": step, "z_steps": zs, "rel_l2": err(tt, o, 1000), "s": dt})
            print(A[-1], flush=True)
    out["z_study"] = A

    C = []
    for rtol in (1e-4, 1e-6, 1e-8, 1e-10):
        tt, zz, o, dt = run(z_steps=400, t_steps=1000, step="ab", rtol=rtol, atol=rtol * 1e-2)
        C.append({"rtol": rtol, "rel_l2": err(tt, o, 1000), "s": dt})
        print(C[-1], flush=True)
    out["rtol_study"] = C

    D = []
    for nt in (250, 500, 1000, 2000):
        tt, zz, o, dt = run(z_steps=400, t_steps=nt, step="ab", rtol=1e-8, atol=1e-10)
        D.append({"t_steps": nt, "rel_l2": err(tt, o, nt), "s": dt})
        print(D[-1], flush=True)
    out["t_study"] = D

    # the reference-problem settings (n_z=200, n_t=1000) and its distance from truth
    tt, zz, o, dt = run(z_steps=200, t_steps=1000, step="ab", rtol=1e-8, atol=1e-10)
    out["reference"] = {"z_steps": 200, "t_steps": 1000, "step": "ab", "rtol": 1e-8,
                        "rel_l2_vs_truth": err(tt, o, 1000), "seconds": dt}
    np.savez("refs/two_level_ref.npz", t=tt, z=zz, omega_zt=o, **{"problem": json.dumps(PROBLEM)})
    json.dump(out, open("refs/convergence.json", "w"), indent=1)
    print("reference", out["reference"])


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "g_scan"
    if which == "convergence":
        convergence()
    if which == "g_scan":
        for g in (0.5, 1.0, 2.0, 5.0):
            t, z, om, dt = run(z_steps=100, problem={"g": g}, t_steps=500)
            e_in, e_out = np.trapezoid(om[0] ** 2, t), np.trapezoid(om[-1] ** 2, t)
            print(f"g={g:4}: peak in {np.abs(om[0]).max():.3f} out {np.abs(om[-1]).max():.3f}"
                  f"  energy out/in {e_out / e_in:.3f}   ({dt:.1f}s)", flush=True)
