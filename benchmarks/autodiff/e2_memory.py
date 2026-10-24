"""E2: peak memory of one reverse-mode gradient against n_z and the number of velocity classes.

Each cell runs in its own process (so the peak resident set is that cell's alone) and records
the resident set before and after, and XLA's compile-time temporary memory when available.

Run:  uv run python benchmarks/autodiff/e2_memory.py --tag e2
"""

import argparse
import json
import resource
import subprocess
import sys
import time

from common import RESULTS, emit

CHECKPOINTS = 32


def rss_mib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20  # bytes on macOS


def child(cfg: dict):
    import jax

    jax.config.update("jax_enable_x64", True)
    from common import Problem, lambda_eit, rb87, time_once

    build = {
        "lambda_eit": lambda: lambda_eit(cfg["n_t"], cfg["n_z"]),
        "rb87": lambda: rb87(cfg["n_t"], cfg["n_z"], cfg["classes_steps"]),
    }[cfg["system"]]
    solver = build()
    opts = {"strategy": cfg["strategy"]}
    if cfg["strategy"] == "s4":
        opts["checkpoints"] = CHECKPOINTS
    prob = Problem(solver, cfg["N"], **opts)
    c = prob.coef0
    before = rss_mib()
    rec = {**cfg, "n_classes": len(solver.thermal_delta_list), "rss_before_mib": before}
    try:
        compiled = jax.jit(jax.value_and_grad(prob.loss)).lower(c).compile()
        ma = compiled.memory_analysis()
        rec["xla_temp_mib"] = (
            getattr(ma, "temp_size_in_bytes", None) and ma.temp_size_in_bytes / 2**20
        )
    except Exception as e:  # memory analysis is not available on every backend
        rec["xla_temp_error"] = repr(e)
    t = time.perf_counter()
    time_once(lambda: prob.value_and_grad(c))
    rec["first_call_s"] = time.perf_counter() - t
    rec["steady_s"] = min(time_once(lambda: prob.value_and_grad(c)) for _ in range(2))
    rec["rss_peak_mib"] = rss_mib()
    print(json.dumps(rec))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tag", default="e2")
    p.add_argument("--quick", action="store_true")
    p.add_argument("--timeout", type=int, default=900)
    p.add_argument("--child", default=None)
    a = p.parse_args()
    if a.child:
        return child(json.loads(a.child))
    cells = []
    n_zs = [20, 40] if a.quick else [50, 100, 200, 400, 800]
    for system, steps in (("lambda_eit", 0), ("rb87", 9)):
        for n_z in n_zs:
            for strat in ("s2", "s3", "s4"):
                cells.append(
                    {
                        "system": system,
                        "n_t": 60 if a.quick else 200,
                        "n_z": n_z,
                        "classes_steps": steps,
                        "strategy": strat,
                        "N": 10,
                        "sweep": "n_z",
                    }
                )
    for steps in [0, 9] if a.quick else [0, 9, 24]:
        for strat in ("s2", "s3", "s4"):
            cells.append(
                {
                    "system": "rb87",
                    "n_t": 60 if a.quick else 200,
                    "n_z": 20 if a.quick else 100,
                    "classes_steps": steps,
                    "strategy": strat,
                    "N": 10,
                    "sweep": "n_v",
                }
            )
    for cell in cells:
        if cell["system"] == "rb87" and cell["strategy"] == "s2" and cell["n_z"] > 400:
            emit(
                {
                    "exp": a.tag,
                    **cell,
                    "skipped": "s2 stores every slice; too large to run here",
                },
                RESULTS / f"{a.tag}.jsonl",
            )
            continue
        cmd = [sys.executable, __file__, "--child", json.dumps(cell)]
        try:
            r = subprocess.run(
                cmd, capture_output=True, text=True, timeout=a.timeout, check=True
            )
            rec = json.loads(r.stdout.strip().splitlines()[-1])
            emit({"exp": a.tag, **rec}, RESULTS / f"{a.tag}.jsonl")
            print(
                cell["system"],
                cell["n_z"],
                cell["classes_steps"],
                cell["strategy"],
                round(rec["rss_peak_mib"]),
                "MiB",
                flush=True,
            )
        except (
            Exception
        ) as e:  # record failures (out of memory, timeout) instead of hiding them
            emit(
                {"exp": a.tag, **cell, "error": repr(e)[:300]},
                RESULTS / f"{a.tag}.jsonl",
            )
            print(cell, "FAILED", repr(e)[:100], flush=True)


if __name__ == "__main__":
    main()
