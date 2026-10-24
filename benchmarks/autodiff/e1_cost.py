"""E1 / E6: the cost of a gradient against the number of parameters.

For each system and parameter count N, times (all warm, compilation excluded):
  forward   one propagation
  reverse   reverse-mode gradient (value_and_grad, strategy s3)
  forward_ad  forward-mode gradient (jacfwd; skipped when it would take too long)
  fd_jax    N + 1 forward evaluations (forward differences) through the compiled propagation
  clerq_solve  one PropagationSolver.solve(); a finite-difference gradient through clerq
              costs (N + 1) of these, reported as an estimate, not measured
Repeats are interleaved across methods so drift shows up in every method alike.

Run:  uv run python benchmarks/autodiff/e1_cost.py --systems two_level lambda_eit rb87 --tag e1
E6 is the same script run against an older source tree (PYTHONPATH=<tree>/src) with --tag e6_before.
"""

import argparse
import time

import numpy as np
from common import RESULTS, Problem, emit, make_solver, summarise, time_once


def make_methods(prob, coef):
    """The callables to time for one problem: the forward pass and the gradient methods."""

    def fwd():
        return prob.ad.propagate({"omega0": prob._field(coef)})

    fd = prob.fd_jax()
    methods = {
        "forward": fwd,
        "reverse": lambda: prob.value_and_grad(coef),
        "fd_jax": lambda: fd(coef),
    }
    return fwd, methods


def forward_mode_method(prob, coef):
    grad = prob.forward_mode_grad()
    return lambda: grad(coef)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--systems", nargs="+", default=["two_level", "lambda_eit", "rb87"])
    p.add_argument("--ns", nargs="+", type=int, default=[2, 5, 10, 50, 100])
    p.add_argument("--repeats", type=int, default=7)
    p.add_argument("--tag", default="e1")
    p.add_argument(
        "--quick", action="store_true", help="tiny grids, 2 repeats (smoke test)"
    )
    p.add_argument("--out", default=None)
    a = p.parse_args()

    import jax

    jax.config.update("jax_enable_x64", True)
    out = a.out or RESULTS / f"{a.tag}.jsonl"
    repeats = 2 if a.quick else a.repeats

    for name in a.systems:
        solver = make_solver(name, a.quick)
        # one clerq solve = the unit of a finite-difference gradient through clerq
        solves = []
        for _ in range(1 if a.quick else 3):
            t = time.perf_counter()
            solver.solve(recalc=True, progress=False, check_stability=False)
            solves.append(time.perf_counter() - t)
        emit(
            {
                "exp": a.tag,
                "system": name,
                "method": "clerq_solve",
                "N": 0,
                **summarise(solves),
            },
            out,
        )

        for n in a.ns:
            prob = Problem(solver, n)
            coef = prob.coef0
            fwd, methods = make_methods(prob, coef)
            skipped = {}
            # time one warm forward pass to size the slow methods
            for _ in range(2):
                time_once(fwd)
            t_fwd = min(time_once(fwd) for _ in range(2))
            if 2.0 * n * t_fwd < 120.0:
                methods["forward_ad"] = forward_mode_method(prob, coef)
            else:
                skipped["forward_ad"] = (
                    f"estimated {2.0 * n * t_fwd:.0f} s per gradient"
                )
            for k, why in skipped.items():
                emit(
                    {"exp": a.tag, "system": name, "method": k, "N": n, "skipped": why},
                    out,
                )
            # warm every method once (compilation), then interleave the repeats
            t_compile = {}
            for k, f in methods.items():
                t = time.perf_counter()
                time_once(f)
                t_compile[k] = time.perf_counter() - t
            slow = sum(t_compile.values()) > 400
            reps = 3 if (slow and not a.quick) else repeats
            times = {k: [] for k in methods}
            for _ in range(reps):
                for k, f in methods.items():
                    times[k].append(time_once(f))
            for k, ts in times.items():
                rec = {
                    "exp": a.tag,
                    "system": name,
                    "method": k,
                    "N": n,
                    "dim": solver.system.num_states,
                    "n_t": len(solver.tlist),
                    "n_z": solver.z_steps,
                    "n_classes": len(solver.thermal_delta_list),
                    "first_call_s": t_compile[k],
                    **summarise(ts),
                }
                emit(rec, out)
            med = {k: float(np.median(v)) for k, v in times.items()}
            print(name, n, {k: round(v, 3) for k, v in med.items()}, flush=True)


if __name__ == "__main__":
    main()
