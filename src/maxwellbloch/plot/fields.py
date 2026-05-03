"""Field plot primitives: field_spacetime, field_envelope, pulse_area."""

from __future__ import annotations

from typing import TYPE_CHECKING

import plotly.graph_objects as go

from maxwellbloch.plot import theme as _theme  # noqa: F401 — registers template
from maxwellbloch.plot.utils import field_label, omega_abs, pulse_area_z

if TYPE_CHECKING:
    from maxwellbloch.mb_solve import MBSolve

_TEMPLATE = "maxwellbloch"


def field_spacetime(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    colorscale: str = "Viridis",
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Heatmap of |Ω(z, t)| for *field_idx*.

    Args:
        mbs: A solved MBSolve instance.
        field_idx: Index of the field to plot.
        colorscale: Plotly colorscale name.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure — call ``.show()`` or ``.write_image()`` on it.
    """
    Omega = omega_abs(mbs, field_idx)
    label = field_label(mbs, field_idx)

    fig = go.Figure(
        data=go.Heatmap(
            z=Omega,
            x=mbs.tlist,
            y=mbs.zlist,
            colorscale=colorscale,
            colorbar=dict(title="|Ω| (γ)"),
        ),
        layout=go.Layout(
            title=f"|Ω(z, t)| — {label}",
            xaxis_title="t (γ⁻¹)",
            yaxis_title="z",
            width=width,
            height=height,
            template=_TEMPLATE,
        ),
    )
    return fig


def field_envelope(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    z_indices: int | list[int] | None = None,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Line plot of |Ω(t)| at one or more fixed z positions.

    Args:
        mbs: A solved MBSolve instance.
        field_idx: Index of the field to plot.
        z_indices: z-step index or list of indices. Defaults to first and last.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    Omega = omega_abs(mbs, field_idx)
    label = field_label(mbs, field_idx)

    if z_indices is None:
        z_indices = [0, -1]
    if isinstance(z_indices, int):
        z_indices = [z_indices]

    fig = go.Figure(
        layout=go.Layout(
            title=f"|Ω(t)| — {label}",
            xaxis_title="t (γ⁻¹)",
            yaxis_title="|Ω| (γ)",
            width=width,
            height=height,
            template=_TEMPLATE,
        )
    )
    for zi in z_indices:
        z_val = mbs.zlist[zi]
        fig.add_trace(
            go.Scatter(x=mbs.tlist, y=Omega[zi], mode="lines", name=f"z = {z_val:.2f}")
        )
    return fig


def pulse_area(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Line plot of pulse area ∫|Ω|dt / π vs z.

    Args:
        mbs: A solved MBSolve instance.
        field_idx: Index of the field to plot.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    area = pulse_area_z(mbs, field_idx)
    label = field_label(mbs, field_idx)

    fig = go.Figure(
        data=go.Scatter(x=mbs.zlist, y=area, mode="lines+markers"),
        layout=go.Layout(
            title=f"Pulse area vs z — {label}",
            xaxis_title="z",
            yaxis_title="area (π)",
            width=width,
            height=height,
            template=_TEMPLATE,
        ),
    )
    return fig
