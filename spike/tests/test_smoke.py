"""T0.1: environment smoke test. diffrax solve + jax.grad in float64 vs analytic."""

import diffrax as dfx
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)


def test_x64_enabled():
    assert jnp.zeros(1).dtype == jnp.float64


def test_grad_matches_analytic():
    # y' = -k y, y(0)=1  ->  y(T) = exp(-k T), dy/dk = -T exp(-k T)
    T = 2.0

    def y_end(k, adjoint=dfx.RecursiveCheckpointAdjoint()):
        sol = dfx.diffeqsolve(
            dfx.ODETerm(lambda t, y, k: -k * y),
            dfx.Tsit5(),
            0.0,
            T,
            0.01,
            jnp.array(1.0),
            args=k,
            stepsize_controller=dfx.PIDController(rtol=1e-10, atol=1e-12),
            adjoint=adjoint,
        )
        return sol.ys[-1] if sol.ys.ndim else sol.ys

    k = 0.7
    g = jax.grad(y_end)(k)
    assert jnp.isclose(g, -T * jnp.exp(-k * T), rtol=1e-8)
    # RecursiveCheckpointAdjoint is reverse-only; forward mode needs ForwardMode()
    fwd = jax.jacfwd(lambda k: y_end(k, dfx.ForwardMode()))(k)
    assert jnp.isclose(fwd, g, rtol=1e-8)
