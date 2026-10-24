"""E8: float32 against float64.

Each problem is run in a fresh process with ``JAX_ENABLE_X64=1`` (the reference) and ``JAX_ENABLE_X64=0`` (float32), both at
the default tolerances for their precision (1e-8 and 1e-10 in float64, 1e-6 and 1e-9 in float32). Reports the error of the
propagated field relative to its peak, of the loss, and of the gradient (relative L2), and the time of one forward pass and
one gradient in each precision.

Run:  uv run python benchmarks/autodiff/e8_precision.py --tag e8
"""

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
from common import RESULTS, emit

CASES = {
    "two_level": ("two_level", {}),
    "lambda_eit": ("lambda_eit", {}),
    "rb87": ("rb87", {}),
    "lambda_eit_long_z": ("lambda_eit", {"n_t": 1000, "n_z": 400}),
    "lambda_eit_fine_t": ("lambda_eit", {"n_t": 4000, "n_z": 100}),
}


def child(case: str, out: str, n_coef: int):
    import common
    import jax

    system, kw = CASES[case]
    prob = common.Problem(getattr(common, system)(**kw), n_coef)
    c = prob.coef0
    field = np.asarray(prob.ad.propagate({"omega0": prob._field(c)}))
    value, grad = prob.value_and_grad(c)
    t_fwd = min(
        common.time_once(lambda: prob.ad.propagate({"omega0": prob._field(c)}))
        for _ in range(3)
    )
    t_grad = min(common.time_once(lambda: prob.value_and_grad(c)) for _ in range(3))
    np.savez(
        out,
        field=field,
        value=float(value),
        grad=np.asarray(grad, dtype=float),
        t_fwd=t_fwd,
        t_grad=t_grad,
        dtype=str(field.dtype),
        x64=bool(jax.config.jax_enable_x64),
    )


def run_child(case: str, x64: bool, n_coef: int) -> dict:
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "r.npz")
        env = {**os.environ, "JAX_ENABLE_X64": "1" if x64 else "0"}
        subprocess.run(
            [
                sys.executable,
                __file__,
                "--child",
                case,
                "--out",
                out,
                "--n",
                str(n_coef),
            ],
            check=True,
            env=env,
        )
        with np.load(out) as z:
            return {k: z[k] for k in z.files}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cases", nargs="+", default=list(CASES))
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--tag", default="e8")
    p.add_argument("--child")
    p.add_argument("--out")
    a = p.parse_args()
    if a.child:
        return child(a.child, a.out, a.n)
    for case in a.cases:
        ref, f32 = run_child(case, True, a.n), run_child(case, False, a.n)
        assert str(ref["dtype"]) == "complex128" and str(f32["dtype"]) == "complex64", (
            ref["dtype"],
            f32["dtype"],
        )
        peak = np.abs(ref["field"]).max()
        rec = {
            "exp": a.tag,
            "case": case,
            "N": a.n,
            "field_error": float(
                np.abs(f32["field"].astype(complex) - ref["field"]).max() / peak
            ),
            "value_error": float(
                abs(float(f32["value"]) - float(ref["value"]))
                / abs(float(ref["value"]))
            ),
            "gradient_error": float(
                np.linalg.norm(f32["grad"] - ref["grad"]) / np.linalg.norm(ref["grad"])
            ),
            "t_fwd_f64": float(ref["t_fwd"]),
            "t_fwd_f32": float(f32["t_fwd"]),
            "t_grad_f64": float(ref["t_grad"]),
            "t_grad_f32": float(f32["t_grad"]),
        }
        emit(rec, RESULTS / f"{a.tag}.jsonl")
        print(
            case,
            {
                k: round(v, 8) if isinstance(v, float) else v
                for k, v in rec.items()
                if k not in ("exp",)
            },
            flush=True,
        )


if __name__ == "__main__":
    main()
