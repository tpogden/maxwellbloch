"""State plot primitives: population, population_spacetime, coherence."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import plotly.graph_objects as go

from maxwellbloch.plot import theme as _theme  # noqa: F401 — registers template

if TYPE_CHECKING:
    from maxwellbloch.mb_solve import MBSolve

_TEMPLATE = "maxwellbloch"


def population(
    mbs: MBSolve,
    state_indices: int | list[int] | None = None,
    *,
    z_idx: int = 0,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Line plot of populations ρ_ii(t) at a fixed z position.

    Args:
        mbs: A solved MBSolve instance.
        state_indices: State index or list of indices. Defaults to all states.
        z_idx: z-step index.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    num_states = mbs.atom.num_states
    if state_indices is None:
        state_indices = list(range(num_states))
    if isinstance(state_indices, int):
        state_indices = [state_indices]

    fig = go.Figure(
        layout=go.Layout(
            title=f"Population vs t (z = {mbs.zlist[z_idx]:.2f})",
            xaxis_title="t (γ⁻¹)",
            yaxis_title="population",
            width=width,
            height=height,
            template=_TEMPLATE,
        )
    )
    for si in state_indices:
        rho_ii = mbs.states_zt[z_idx, :, si, si].real
        fig.add_trace(go.Scatter(x=mbs.tlist, y=rho_ii, mode="lines", name=f"|{si}⟩"))
    fig.update_layout(
        xaxis_minallowed=float(mbs.tlist[0]),
        xaxis_maxallowed=float(mbs.tlist[-1]),
    )
    return fig


def population_spacetime(
    mbs: MBSolve,
    state_idx: int = 0,
    *,
    colorscale: str = "Viridis",
    plot_type: Literal["heatmap", "contour"] = "heatmap",
    zmin: float | None = None,
    zmax: float | None = None,
    show_z_bounds: tuple[float, float] | None = None,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Heatmap or contour plot of population ρ_ii(z, t) for a single state.

    Args:
        mbs: A solved MBSolve instance.
        state_idx: State index.
        colorscale: Plotly colorscale name.
        plot_type: ``"heatmap"`` (default, smoothed) or ``"contour"``.
        zmin: Colorbar minimum. Defaults to auto (data minimum).
        zmax: Colorbar maximum. Defaults to auto (data maximum). Pass
            ``zmin=0, zmax=1`` to compare populations on a common scale.
        show_z_bounds: If given, draw dotted grey lines at these two z values
            to mark the start and end of the medium, e.g. ``(0.0, 1.0)``.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    rho_ii = mbs.states_zt[:, :, state_idx, state_idx].real

    _common = dict(z=rho_ii, x=mbs.tlist, y=mbs.zlist, colorscale=colorscale,
                   colorbar=dict(title=f"ρ_{state_idx}{state_idx}"),
                   zmin=zmin, zmax=zmax)
    if plot_type == "contour":
        trace = go.Contour(**_common, line_width=0.5, line_smoothing=0.85)
    else:
        trace = go.Heatmap(**_common, zsmooth="best")

    fig = go.Figure(
        data=trace,
        layout=go.Layout(
            title=f"Population ρ_{state_idx}{state_idx}(z, t)",
            xaxis_title="t (γ⁻¹)",
            yaxis_title="z",
            width=width,
            height=height,
            template=_TEMPLATE,
        ),
    )
    fig.update_layout(
        xaxis_minallowed=float(mbs.tlist[0]),
        xaxis_maxallowed=float(mbs.tlist[-1]),
        yaxis_minallowed=float(mbs.zlist[0]),
        yaxis_maxallowed=float(mbs.zlist[-1]),
    )
    if show_z_bounds is not None:
        for z in show_z_bounds:
            fig.add_hline(y=float(z), line=dict(color="grey", dash="dot", width=1))
    return fig


def coherence(
    mbs: MBSolve,
    i: int,
    j: int,
    *,
    z_idx: int = 0,
    component: str = "abs",
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Line plot of coherence ρ_ij(t) at a fixed z position.

    Args:
        mbs: A solved MBSolve instance.
        i: Row index of the density matrix element.
        j: Column index of the density matrix element.
        z_idx: z-step index.
        component: One of ``"abs"``, ``"real"``, or ``"imag"``.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    rho_ij = mbs.states_zt[z_idx, :, i, j]
    if component == "abs":
        y = np.abs(rho_ij)
        ylabel = f"|ρ_{i}{j}|"
    elif component == "real":
        y = rho_ij.real
        ylabel = f"Re(ρ_{i}{j})"
    elif component == "imag":
        y = rho_ij.imag
        ylabel = f"Im(ρ_{i}{j})"
    else:
        raise ValueError(
            f"component must be 'abs', 'real', or 'imag', got {component!r}"
        )

    fig = go.Figure(
        data=go.Scatter(x=mbs.tlist, y=y, mode="lines"),
        layout=go.Layout(
            title=f"{ylabel}(t) at z = {mbs.zlist[z_idx]:.2f}",
            xaxis_title="t (γ⁻¹)",
            yaxis_title=ylabel,
            width=width,
            height=height,
            template=_TEMPLATE,
        ),
    )
    fig.update_layout(
        xaxis_minallowed=float(mbs.tlist[0]),
        xaxis_maxallowed=float(mbs.tlist[-1]),
    )
    return fig
