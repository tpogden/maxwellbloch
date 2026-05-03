"""Spectral plot primitives: spectrum, spectrum_overlay."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import plotly.graph_objects as go

from maxwellbloch import spectral
from maxwellbloch.plot import theme as _theme  # noqa: F401 — registers template

if TYPE_CHECKING:
    from maxwellbloch.mb_solve import MBSolve

_TEMPLATE = "maxwellbloch"

# Candidate tick magnitudes for arcsinh axis (symmetric about 0).
# Ticks ≤ the displayed f_max are shown; no hardcoded upper bound.
_ARCSINH_TICK_CANDIDATES = [1, 2, 5, 10, 20, 50, 100, 200]


def _apply_freq_options(
    freqs: np.ndarray,
    arrays: list[np.ndarray],
    freq_range: float | None,
    freq_scale: str,
    arcsinh_scale: float,
) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
    """Clip and/or transform the frequency axis and matching data arrays.

    Returns ``(x_plot, clipped_arrays, freqs_clipped)`` where ``x_plot`` is
    the array of x-values to pass to Plotly traces.
    """
    if freq_range is not None:
        mask = np.abs(freqs) <= freq_range
        freqs = freqs[mask]
        arrays = [a[mask] for a in arrays]

    if freq_scale == "arcsinh":
        x_plot = np.arcsinh(freqs / arcsinh_scale)
    else:
        x_plot = freqs

    return x_plot, arrays, freqs


def _arcsinh_tick_layout(freqs: np.ndarray, arcsinh_scale: float) -> dict:
    """Return an xaxis dict with human-readable tick labels for an arcsinh axis."""
    f_max = np.max(np.abs(freqs))
    tick_vals = sorted(
        {0}
        | {v for v in _ARCSINH_TICK_CANDIDATES if v <= f_max}
        | {-v for v in _ARCSINH_TICK_CANDIDATES if v <= f_max}
    )
    return dict(
        tickmode="array",
        tickvals=np.arcsinh(np.array(tick_vals) / arcsinh_scale).tolist(),
        ticktext=[str(v) for v in tick_vals],
    )


def spectrum(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    z_idx: int = -1,
    show_dispersion: bool = False,
    freq_range: float | None = None,
    freq_scale: str = "linear",
    arcsinh_scale: float = 1.0,
    window: str | None = None,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Absorption (and optionally dispersion) spectrum via Fourier transform.

    Uses the Fourier transform of the time-domain field to extract the
    frequency-domain absorption at z position *z_idx*.

    Args:
        mbs: A solved MBSolve instance. The probe field should be a weak
            pulse (so the linear-response approximation holds).
        field_idx: Index of the probe field.
        z_idx: z-step index at which to compute absorption.
        show_dispersion: If True, also plot the dispersion on a secondary axis.
        freq_range: If given, only frequencies with |f| ≤ freq_range are shown.
            Clips the underlying data arrays, not just the axis view.
        freq_scale: ``"linear"`` (default) or ``"arcsinh"``. With ``"arcsinh"``
            the x-axis is compressed nonlinearly — linear near zero, log-like at
            large |f| — so narrow and broad features coexist in one plot.
        arcsinh_scale: Scale factor for the arcsinh transform (only used when
            *freq_scale* is ``"arcsinh"``). Larger values give a wider linear
            region near zero.
        window: Optional SciPy window name (e.g. ``"hann"``, ``"blackman"``)
            applied to the time-domain field before the FFT. Reduces spectral
            leakage from step-function pulses. Default ``None`` (no windowing).
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    freqs = spectral.freq_list(mbs)
    arrays = [spectral.absorption(mbs, field_idx, z_idx=z_idx, window=window)]
    if show_dispersion:
        arrays.append(spectral.dispersion(mbs, field_idx, z_idx=z_idx, window=window))

    label = mbs.atom.fields[field_idx].label or f"field {field_idx}"

    x_plot, arrays, freqs_clipped = _apply_freq_options(
        freqs, arrays, freq_range, freq_scale, arcsinh_scale
    )
    absorp = arrays[0]

    fig = go.Figure(
        layout=go.Layout(
            title=f"Absorption spectrum — {label}",
            xaxis_title="frequency (γ)",
            yaxis_title="absorption (a.u.)",
            width=width,
            height=height,
            template=_TEMPLATE,
        )
    )
    fig.add_trace(go.Scatter(x=x_plot, y=absorp, mode="lines", name="absorption"))

    if show_dispersion:
        fig.add_trace(
            go.Scatter(
                x=x_plot,
                y=arrays[1],
                mode="lines",
                name="dispersion",
                yaxis="y2",
            )
        )
        fig.update_layout(
            yaxis2=dict(
                title="dispersion (a.u.)",
                overlaying="y",
                side="right",
                showgrid=False,
            )
        )

    if freq_scale == "arcsinh":
        fig.update_layout(xaxis=_arcsinh_tick_layout(freqs_clipped, arcsinh_scale))

    return fig


def spectrum_overlay(
    mbs_list: list[MBSolve],
    field_idx: int = 0,
    *,
    z_idx: int = -1,
    labels: list[str] | None = None,
    freq_range: float | None = None,
    freq_scale: str = "linear",
    arcsinh_scale: float = 1.0,
    window: str | None = None,
    width: int = 900,
    height: int = 400,
) -> go.Figure:
    """Overlay absorption spectra from multiple MBSolve results.

    Args:
        mbs_list: List of solved MBSolve instances.
        field_idx: Field index to plot in each.
        z_idx: z-step index for absorption.
        labels: Trace labels, one per element of mbs_list.
        freq_range: If given, only frequencies with |f| ≤ freq_range are shown.
            Clips the underlying data arrays, not just the axis view.
        freq_scale: ``"linear"`` (default) or ``"arcsinh"``. With ``"arcsinh"``
            the x-axis is compressed nonlinearly — linear near zero, log-like at
            large |f|.
        arcsinh_scale: Scale factor for the arcsinh transform (only used when
            *freq_scale* is ``"arcsinh"``).
        window: Optional SciPy window name (e.g. ``"hann"``, ``"blackman"``)
            applied to the time-domain field before the FFT. Default ``None``.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    if labels is None:
        labels = [f"run {i}" for i in range(len(mbs_list))]

    fig = go.Figure(
        layout=go.Layout(
            title="Absorption spectra",
            xaxis_title="frequency (γ)",
            yaxis_title="absorption (a.u.)",
            width=width,
            height=height,
            template=_TEMPLATE,
        )
    )
    freqs_clipped_last = None
    for mbs, lbl in zip(mbs_list, labels, strict=False):
        freqs = spectral.freq_list(mbs)
        absorp = spectral.absorption(mbs, field_idx, z_idx=z_idx, window=window)
        x_plot, (absorp,), freqs_clipped = _apply_freq_options(
            freqs, [absorp], freq_range, freq_scale, arcsinh_scale
        )
        freqs_clipped_last = freqs_clipped
        fig.add_trace(go.Scatter(x=x_plot, y=absorp, mode="lines", name=lbl))

    if freq_scale == "arcsinh" and freqs_clipped_last is not None:
        fig.update_layout(xaxis=_arcsinh_tick_layout(freqs_clipped_last, arcsinh_scale))

    return fig
