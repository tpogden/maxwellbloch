"""T0.5 + Phase-2-style check on the toy chain: AD gradients vs FD and vs jacfwd.

Fixed-step Tsit5 (ConstantStepSize) so the discretised map is smooth and FD is meaningful.
Gradient is w.r.t. the input field samples, through n_z interpolated-drive ODE solves.
"""

import sys

import jax
import jax.numpy as jnp
import numpy as np
import pytest

sys.path.insert(0, "bench")
from probe_memory import make_problem  # noqa: E402

from mbjax.fd import check_gradient, rel_err  # noqa: E402

jax.config.update("jax_enable_x64", True)

N_Z, N_T, DIM = 6, 60, 4
IDX = [0, 7, 19, 30, 31, 44, 59]


@pytest.mark.parametrize("strategy", ["direct", "s2", "s3", "s5"])
def test_reverse_vs_fd(strategy):
    loss, x0 = make_problem(N_Z, N_T, DIM, strategy)
    r = check_gradient(loss, x0, IDX)
    assert min(r["rel_err"].values()) < 1e-5, r["rel_err"]


def test_reverse_matches_forward_mode():
    loss_r, x0 = make_problem(N_Z, N_T, DIM, "s3")
    loss_f, _ = make_problem(N_Z, N_T, DIM, "fwd")
    g_r = jax.jit(jax.grad(loss_r))(x0)
    g_f = jax.jit(jax.jacfwd(loss_f))(x0)
    assert rel_err(g_r, g_f) < 1e-10


def test_strategies_agree_exactly():
    gs = {}
    for s in ("direct", "s2", "s3", "s5"):
        loss, x0 = make_problem(N_Z, N_T, DIM, s)
        gs[s] = np.asarray(jax.jit(jax.grad(loss))(x0))
    for s in ("s2", "s3", "s5"):
        assert rel_err(gs[s], gs["direct"]) < 1e-10
