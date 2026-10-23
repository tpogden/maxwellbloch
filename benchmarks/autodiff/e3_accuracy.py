"""E3: accuracy of the gradients and of the propagated field.

Per system (moderate grids, so fixed-step references are affordable), with N coefficients:
  field_error      max |propagate - solve()| / peak, at the problem's own parameters
  ad_vs_fd_fixed   reverse AD vs central finite differences (best of four steps), both on a fixed-step solve
                   (fixed steps, not adaptive: adaptive finite differences are noise limited)
  fwd_vs_rev       forward-mode vs reverse-mode gradient (relative L2 difference)
  adaptive_vs_fixed  reverse AD with the adaptive solver vs reverse AD on the fixed-step solve
"""

import argparse

import numpy as np
from common import RESULTS, Problem, emit, two_level, lambda_eit, rb87

# system: (builder, (n_t, n_z), N coefficients, fixed steps per slice); the step counts are not multiples of
# the output intervals, and the Rb87 case is smaller because 11 classes on a fixed grid are slow
GRIDS = {
    "two_level": (two_level, (200, 40), 10, 2003),
    "lambda_eit": (lambda_eit, (100, 20), 10, 2003),
    "rb87": (rb87, (100, 20), 6, 1009),
}


def rel(a, b):
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)) / max(np.linalg.norm(np.asarray(b)), 1e-300))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--systems", nargs="+", default=list(GRIDS))
    p.add_argument("--tag", default="e3")
    p.add_argument("--quick", action="store_true")
    a = p.parse_args()

    import jax

    jax.config.update("jax_enable_x64", True)
    out = RESULTS / f"{a.tag}.jsonl"
    for name in a.systems:
        build, (n_t, n_z), n_coef, n_fixed = GRIDS[name]
        solver = build(n_t // 2, n_z // 2) if a.quick else build(n_t, n_z)
        solver.solve(recalc=True, progress=False, check_stability=False)
        ref = np.asarray(solver.rabi_freqs_zt)
        adaptive = Problem(solver, n_coef, rtol=1e-8, atol=1e-10)
        got = np.asarray(adaptive.ad.propagate())
        field_error = float(np.abs(got - ref).max() / np.abs(ref).max())

        fixed = Problem(solver, n_coef, n_steps=n_fixed)
        c0 = fixed.coef0
        _, g_fixed = fixed.value_and_grad(c0)
        # central differences at several steps: truncation error falls and rounding error rises
        # as h shrinks, so report the whole sequence and its best value
        by_h = {}
        for h in (1e-5, 1e-6, 1e-7, 1e-8):
            fd = np.array([(float(fixed.value(c0 + h * np.eye(n_coef)[i])) - float(fixed.value(c0 - h * np.eye(n_coef)[i]))) / (2 * h)
                           for i in range(n_coef)])
            by_h[f"{h:g}"] = rel(g_fixed, fd)
        fd_norm = float(np.linalg.norm(np.asarray(g_fixed)))
        _, g_adapt = adaptive.value_and_grad(c0)
        _, g_fwd = fixed.forward_mode_grad()(c0)
        emit({"exp": a.tag, "system": name, "N": n_coef, "n_steps_fixed": n_fixed, "n_t": n_t, "n_z": n_z, "n_classes": len(solver.thermal_delta_list),
              "field_error_vs_solve": field_error, "ad_vs_fd_fixed": min(by_h.values()), "ad_vs_fd_by_step": by_h,
              "fwd_vs_rev": rel(g_fwd, g_fixed), "adaptive_vs_fixed": rel(g_adapt, g_fixed),
              "grad_norm": fd_norm}, out)
        print(name, "done", flush=True)


if __name__ == "__main__":
    main()
