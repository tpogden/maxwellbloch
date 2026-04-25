"""MaxwellBloch Plotly theme template."""

import plotly.graph_objects as go
import plotly.io as pio

_CATEGORICAL_COLORS = [
    "#4477AA",  # blue
    "#EE6677",  # rose
    "#228833",  # green
    "#CCBB44",  # yellow
    "#66CCEE",  # cyan
    "#AA3377",  # purple
    "#BBBBBB",  # grey
]

_TEMPLATE = go.layout.Template(
    layout=go.Layout(
        font=dict(family="sans-serif", size=13),
        paper_bgcolor="white",
        plot_bgcolor="white",
        colorway=_CATEGORICAL_COLORS,
        xaxis=dict(
            showgrid=True,
            gridcolor="#e0e0e0",
            zeroline=True,
            zerolinecolor="#aaaaaa",
            zerolinewidth=1,
            ticks="inside",
            showline=True,
            linecolor="#888888",
            mirror=False,
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor="#e0e0e0",
            zeroline=True,
            zerolinecolor="#aaaaaa",
            zerolinewidth=1,
            ticks="inside",
            showline=True,
            linecolor="#888888",
            mirror=False,
        ),
        colorscale=dict(sequential="Viridis", diverging="RdBu"),
        legend=dict(
            bgcolor="rgba(255,255,255,0.8)", bordercolor="#cccccc", borderwidth=1
        ),
        margin=dict(l=60, r=20, t=50, b=60),
    )
)

pio.templates["maxwellbloch"] = _TEMPLATE
