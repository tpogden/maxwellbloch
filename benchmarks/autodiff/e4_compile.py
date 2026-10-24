"""E4: compile time and first-call vs steady-state time (each system in a fresh process).

Run:  uv run python benchmarks/autodiff/e4_compile.py --tag e4
"""

import argparse
import subprocess
import sys
import time

from common import RESULTS, SYSTEMS, Problem, emit, make_solver, time_once


def child(name: str, n: int, quick: bool, tag: str):
    import jax

    jax.config.update("jax_enable_x64", True)
    prob = Problem(make_solver(name, quick), n)
    c = prob.coef0
    rec = {"exp": tag, "system": name, "N": n}
    for label, f in (
        ("forward", lambda: prob.ad.propagate({"omega0": prob._field(c)})),
        ("reverse", lambda: prob.value_and_grad(c)),
    ):
        first = time_once(f)
        steady = min(time_once(f) for _ in range(3))
        rec[f"{label}_first_s"], rec[f"{label}_steady_s"] = first, steady
    emit(rec, RESULTS / f"{tag}.jsonl")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--systems", nargs="+", default=list(SYSTEMS))
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--tag", default="e4")
    p.add_argument("--quick", action="store_true")
    p.add_argument("--child", default=None)
    a = p.parse_args()
    if a.child:
        return child(a.child, a.n, a.quick, a.tag)
    for name in a.systems:
        cmd = [
            sys.executable,
            __file__,
            "--child",
            name,
            "--n",
            str(a.n),
            "--tag",
            a.tag,
        ] + (["--quick"] if a.quick else [])
        t = time.perf_counter()
        subprocess.run(cmd, check=True)
        print(
            name,
            "compile benchmark done in",
            round(time.perf_counter() - t, 1),
            "s",
            flush=True,
        )


if __name__ == "__main__":
    main()
