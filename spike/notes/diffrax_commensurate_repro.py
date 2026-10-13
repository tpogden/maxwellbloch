"""Minimal repro: RecursiveCheckpointAdjoint gives wrong reverse-mode gradients when
SaveAt(ts) times coincide with ConstantStepSize step boundaries (3.5e-2 here, eager and jit;
in a Lindblad slice the same failure appeared only under jit, 2.8x error).

Observed with jax/jaxlib 0.11.2, diffrax 0.7.2, equinox 0.13.8, CPU, float64.
Run: uv run python notes/diffrax_commensurate_repro.py
Expected: every row's AD matches the central difference. Actual: the n=1500
RecursiveCheckpointAdjoint rows are wrong (forward values are fine); n=1501 is fine, and
DirectAdjoint / ForwardMode are fine at both.
"""
import diffrax as dfx
import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

ts = jnp.linspace(0.0, 6.0, 61)  # output spacing 0.1


def make_loss(adjoint, n_steps):
    def f(t, y, th):
        return -th[0] * y + th[1] * jnp.cos(3.0 * t)

    def loss(th):
        sol = dfx.diffeqsolve(
            dfx.ODETerm(f), dfx.Tsit5(), 0.0, 6.0, 6.0 / n_steps, jnp.array(1.0), args=th,
            saveat=dfx.SaveAt(ts=ts), stepsize_controller=dfx.ConstantStepSize(),
            adjoint=adjoint, max_steps=n_steps + 2)
        return jnp.sum(jnp.arange(61.0) * sol.ys)

    return loss


th = jnp.array([1.1, 0.4])
for n in (1500, 1501):
    L = make_loss(dfx.RecursiveCheckpointAdjoint(), n)
    h = 1e-4
    fd = np.array([(L(th + h * e) - L(th - h * e)) / (2 * h) for e in jnp.eye(2)])
    for name, adj, mode in [("RecursiveCheckpointAdjoint", dfx.RecursiveCheckpointAdjoint(), "rev"),
                            ("DirectAdjoint", dfx.DirectAdjoint(), "rev"),
                            ("ForwardMode", dfx.ForwardMode(), "fwd")]:
        for jit in (False, True):
            g = jax.jacfwd(make_loss(adj, n)) if mode == "fwd" else jax.grad(make_loss(adj, n))
            g = np.asarray((jax.jit(g) if jit else g)(th))
            err = np.linalg.norm(g - fd) / np.linalg.norm(fd)
            print(f"n={n} steps/interval={n/60:8.3f} {name:27} jit={jit!s:5} rel err vs FD {err:.1e}")
