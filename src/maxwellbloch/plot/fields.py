"""Field plot primitives: field_spacetime, field_envelope, pulse_area."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import plotly.graph_objects as go

from maxwellbloch.plot import theme as _theme  # noqa: F401 — registers template
from maxwellbloch.plot.utils import field_label, omega_abs, pulse_area_z

if TYPE_CHECKING:
    import numpy as np

    from maxwellbloch.field import Field
    from maxwellbloch.mb_solve import MBSolve

_TEMPLATE = "maxwellbloch"

_DEFAULT_COLORSCALES = ["Blues", "Greens", "Oranges", "Reds"]

# (line color, semitransparent fill) — midpoints of each colorscale above.
_DEFAULT_FIELD_COLORS = [
    ("#2171B5", "rgba(33, 113, 181, 0.2)"),
    ("#238B45", "rgba(35, 139, 69, 0.2)"),
    ("#D94801", "rgba(217, 72, 1, 0.2)"),
    ("#CB181D", "rgba(203, 24, 29, 0.2)"),
]


def field_profile(
    field: "Field",
    tlist: "np.ndarray",
    *,
    field_idx: int = 0,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Line plot of the input Rabi frequency profile |Ω(t)| for a single field.

    Args:
        field: A Field instance with rabi_freq_t_func and rabi_freq_t_args set.
        tlist: 1-D array of time points to evaluate the profile over.
        field_idx: Used only to pick a colour from the default palette.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    import numpy as np

    profile = np.abs(
        field.rabi_freq * field.rabi_freq_t_func(tlist, args=field.rabi_freq_t_args)
    )
    line_color, fill_color = _DEFAULT_FIELD_COLORS[
        field_idx % len(_DEFAULT_FIELD_COLORS)
    ]

    fig = go.Figure(
        data=go.Scatter(
            x=tlist,
            y=profile,
            mode="lines",
            line=dict(color=line_color),
            fill="tozeroy",
            fillcolor=fill_color,
            showlegend=False,
        ),
        layout=go.Layout(
            title="|Ω(t)| — input profile",
            xaxis_title="t",
            yaxis_title="|Ω|",
            width=width,
            height=height,
            template=_TEMPLATE,
        ),
    )
    return fig


def field_spacetime(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    colorscale: str | None = None,
    plot_type: Literal["heatmap", "contour"] = "heatmap",
    show_z_bounds: tuple[float, float] | None = None,
    speed_of_light: float | None = None,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Heatmap or contour plot of |Ω(z, t)| for *field_idx*.

    Args:
        mbs: A solved MBSolve instance.
        field_idx: Index of the field to plot.
        colorscale: Plotly colorscale name. Defaults to Blues for field 0,
            Greens for field 1 (cycles through _DEFAULT_COLORSCALES).
        plot_type: ``"heatmap"`` (default, smoothed) or ``"contour"``.
        show_z_bounds: If given, draw dotted grey lines at these two z values
            to mark the start and end of the medium, e.g. ``(0.0, 1.0)``.
        speed_of_light: If given, apply the fixed-frame transform using
            ``fixed.t_list`` and ``fixed.rabi_freq`` at this speed.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure — call ``.show()`` or ``.write_image()`` on it.
    """
    if colorscale is None:
        colorscale = _DEFAULT_COLORSCALES[field_idx % len(_DEFAULT_COLORSCALES)]
    label = field_label(mbs, field_idx)

    if speed_of_light is not None:
        from maxwellbloch import fixed

        Omega = fixed.rabi_freq(mbs, field_idx, speed_of_light, part="abs")
        tlist = fixed.t_list(mbs, speed_of_light)
        x_title = "t (fixed frame)"
        title_suffix = " — fixed frame"
    else:
        Omega = omega_abs(mbs, field_idx)
        tlist = mbs.tlist
        x_title = "t (γ⁻¹)"
        title_suffix = ""

    _common = dict(
        z=Omega,
        x=tlist,
        y=mbs.zlist,
        colorscale=colorscale,
        colorbar=dict(title="|Ω| (γ)"),
    )
    if plot_type == "contour":
        trace = go.Contour(**_common, line_width=0.5, line_smoothing=0.85)
    else:
        trace = go.Heatmap(**_common, zsmooth="best")

    fig = go.Figure(
        data=trace,
        layout=go.Layout(
            title=f"|Ω(z, t)| — {label}{title_suffix}",
            xaxis_title=x_title,
            yaxis_title="z",
            width=width,
            height=height,
            template=_TEMPLATE,
        ),
    )
    fig.update_layout(
        xaxis_minallowed=float(tlist[0]),
        xaxis_maxallowed=float(tlist[-1]),
        yaxis_minallowed=float(mbs.zlist[0]),
        yaxis_maxallowed=float(mbs.zlist[-1]),
    )
    if show_z_bounds is not None:
        for z in show_z_bounds:
            fig.add_hline(y=float(z), line=dict(color="grey", dash="dot", width=1))
    return fig


def field_envelope(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    z_indices: int | list[int] | None = None,
    speed_of_light: float | None = None,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Line plot of |Ω(t)| at one or more fixed z positions.

    Args:
        mbs: A solved MBSolve instance.
        field_idx: Index of the field to plot.
        z_indices: z-step index or list of indices. Defaults to first and last.
        speed_of_light: If given, apply the fixed-frame transform using
            ``fixed.t_list`` and ``fixed.rabi_freq`` at this speed.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    label = field_label(mbs, field_idx)
    line_color, fill_color = _DEFAULT_FIELD_COLORS[
        field_idx % len(_DEFAULT_FIELD_COLORS)
    ]

    if speed_of_light is not None:
        from maxwellbloch import fixed

        Omega = fixed.rabi_freq(mbs, field_idx, speed_of_light, part="abs")
        tlist = fixed.t_list(mbs, speed_of_light)
        x_title = "t (fixed frame)"
        title_suffix = " — fixed frame"
    else:
        Omega = omega_abs(mbs, field_idx)
        tlist = mbs.tlist
        x_title = "t (γ⁻¹)"
        title_suffix = ""

    if z_indices is None:
        z_indices = [0, -1]
    if isinstance(z_indices, int):
        z_indices = [z_indices]

    fig = go.Figure(
        layout=go.Layout(
            title=f"|Ω(t)| — {label}{title_suffix}",
            xaxis_title=x_title,
            yaxis_title="|Ω| (γ)",
            width=width,
            height=height,
            template=_TEMPLATE,
        )
    )
    for zi in z_indices:
        z_val = mbs.zlist[zi]
        fig.add_trace(
            go.Scatter(
                x=tlist,
                y=Omega[zi],
                mode="lines",
                name=f"z = {z_val:.2f}",
                line=dict(color=line_color),
                fill="tozeroy",
                fillcolor=fill_color,
            )
        )
    fig.update_layout(
        xaxis_minallowed=float(tlist[0]),
        xaxis_maxallowed=float(tlist[-1]),
    )
    return fig


def field_z_profile_anim(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    speed_of_light: float | None = None,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Animated spatial profile |Ω(z)| at a selected time t.

    Displays the field envelope along z at a single time step. Use the slider
    to scrub through time, or press Play to animate the pulse propagation.

    Args:
        mbs: A solved MBSolve instance.
        field_idx: Index of the field to plot.
        speed_of_light: If given, apply the fixed-frame transform using
            ``fixed.t_list`` and ``fixed.rabi_freq`` at this speed.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure with animation frames, time slider, and play/pause buttons.
    """
    if speed_of_light is not None:
        from maxwellbloch import fixed

        Omega = fixed.rabi_freq(mbs, field_idx, speed_of_light, part="abs")
        tlist = fixed.t_list(mbs, speed_of_light)
    else:
        Omega = omega_abs(mbs, field_idx)
        tlist = mbs.tlist

    label = field_label(mbs, field_idx)
    line_color, fill_color = _DEFAULT_FIELD_COLORS[
        field_idx % len(_DEFAULT_FIELD_COLORS)
    ]
    y_max = float(Omega.max())

    frames = [
        go.Frame(data=[go.Scatter(y=Omega[:, t_idx])], name=str(t_idx))
        for t_idx in range(len(tlist))
    ]

    slider_steps = [
        {
            "args": [
                [str(t_idx)],
                {
                    "frame": {"duration": 0, "redraw": False},
                    "mode": "immediate",
                    "transition": {"duration": 0},
                },
            ],
            "label": f"{t_val:.2f}",
            "method": "animate",
        }
        for t_idx, t_val in enumerate(tlist)
    ]

    fig = go.Figure(
        data=[
            go.Scatter(
                x=mbs.zlist,
                y=Omega[:, 0],
                mode="lines",
                line=dict(color=line_color),
                fill="tozeroy",
                fillcolor=fill_color,
                showlegend=False,
            )
        ],
        frames=frames,
        layout=go.Layout(
            title=f"|Ω(z, t)| — {label}",
            xaxis_title="z",
            yaxis_title="|Ω| (γ)",
            yaxis_range=[0, y_max * 1.05],
            width=width,
            height=height,
            template=_TEMPLATE,
            updatemenus=[
                {
                    "buttons": [
                        {
                            "args": [
                                None,
                                {
                                    "frame": {"duration": 50, "redraw": False},
                                    "fromcurrent": True,
                                    "transition": {"duration": 0},
                                },
                            ],
                            "label": "▶ Play",
                            "method": "animate",
                        },
                        {
                            "args": [
                                [None],
                                {
                                    "frame": {"duration": 0, "redraw": False},
                                    "mode": "immediate",
                                    "transition": {"duration": 0},
                                },
                            ],
                            "label": "⏸ Pause",
                            "method": "animate",
                        },
                    ],
                    "direction": "left",
                    "pad": {"r": 10, "t": 87},
                    "showactive": False,
                    "type": "buttons",
                    "x": 0.1,
                    "xanchor": "right",
                    "y": 0,
                    "yanchor": "top",
                }
            ],
            sliders=[
                {
                    "active": 0,
                    "yanchor": "top",
                    "xanchor": "left",
                    "currentvalue": {
                        "prefix": "t = ",
                        "visible": True,
                        "xanchor": "right",
                    },
                    "transition": {"duration": 0},
                    "pad": {"b": 10, "t": 50},
                    "len": 0.9,
                    "x": 0.1,
                    "y": 0,
                    "steps": slider_steps,
                }
            ],
        ),
    )
    fig.update_layout(
        xaxis_minallowed=float(mbs.zlist[0]),
        xaxis_maxallowed=float(mbs.zlist[-1]),
        yaxis_minallowed=0.0,
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

    line_color, fill_color = _DEFAULT_FIELD_COLORS[
        2
    ]  # orange — distinct from field blue

    fig = go.Figure(
        data=go.Scatter(
            x=mbs.zlist,
            y=area,
            mode="lines",
            fill="tozeroy",
            line=dict(color=line_color),
            fillcolor=fill_color,
        ),
        layout=go.Layout(
            title=f"Pulse area vs z — {label}",
            xaxis_title="z",
            yaxis_title="area (π)",
            width=width,
            height=height,
            template=_TEMPLATE,
        ),
    )
    fig.update_layout(
        xaxis_minallowed=float(mbs.zlist[0]),
        xaxis_maxallowed=float(mbs.zlist[-1]),
    )
    return fig
