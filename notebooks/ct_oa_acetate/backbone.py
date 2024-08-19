import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from glob import glob
import pandas as pd
from fitting import *

colors = {
    "ct": "#7570B3",
    "oa": "#D95F02",
}

w = 400
h = 250


def style_plot(fig, marker_size=3):
    """Style function for figures setting fot size and true black color."""
    for d in fig["data"]:
        d["marker"]["size"] = marker_size
        d["line"]["width"] = marker_size
    # Font size
    j = 10
    fig.update_layout(font={"size": j, "color": "black"})
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = j
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = j
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = j
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"
    fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=j, color="black")))
    fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(size=j, color="black")))
    fig.update_layout(
        margin=dict(l=60, r=20, t=20, b=20),
    )
    fig.update_yaxes(title_standoff=10)
    fig.update_xaxes(title_standoff=10)
    return fig


def get_ct_oa():
    meta, raw = get_dfs("./")
    # Filter data
    mask = (meta["species"] == "Comamonas testosteroni") & (meta["cs_conc"] == 15)
    columns = list(set(meta[mask]["linegroup"]))
    ct = raw[["time"] + columns].dropna()
    ct = mask_df(ct, 0, 30)

    meta, raw = get_dfs("./")
    # Filter data
    mask = (meta["species"] == "Ochrobactrum anthropi") & (meta["cs_conc"] == 15)
    columns = list(set(meta[mask]["linegroup"]))
    oa = raw[["time"] + columns].dropna()
    oa = mask_df(oa, 14, 22)

    return ct, oa


def fig1():
    ct, oa = get_ct_oa()
    fig = make_subplots(
        rows=1,
        cols=2,
        vertical_spacing=0.2,
        horizontal_spacing=0.12,
        row_heights=[0.3],
        subplot_titles=["Ct", "Oa"],
    )
    for c in ct.columns[1:]:
        fig.add_trace(
            go.Scatter(
                x=ct["time"],
                y=ct[c],
                marker=dict(color=colors["ct"]),
                mode="lines",
                showlegend=False,
            ),
            row=1,
            col=1,
        )
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")

    for c in oa.columns[1:]:
        fig.add_trace(
            go.Scatter(
                x=oa["time"],
                y=oa[c],
                marker=dict(color=colors["oa"]),
                mode="lines",
                showlegend=False,
            ),
            row=1,
            col=2,
        )
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig = style_plot(fig)
    fig.show()


def fig2():
    ct, oa = get_ct_oa()
    t1, t2 = 0, 5
    vct = fit_max_growth_rate(ct, t1, t2, plot=False)
    ct = ct[(ct["time"] >= t1) & (ct["time"] <= t2)]
    n0 = get_n0(ct)
    fit = [n0 * np.exp(vct * t) for t in ct["time"]]
    fig = make_subplots(
        rows=1,
        cols=2,
        vertical_spacing=0.2,
        horizontal_spacing=0.12,
        row_heights=[0.3],
        subplot_titles=["Ct", "Oa"],
    )
    for c in ct.columns[1:]:
        fig.add_trace(
            go.Scatter(
                x=ct["time"],
                y=ct[c],
                marker=dict(color=colors["ct"]),
                mode="lines",
                showlegend=False,
            ),
            row=1,
            col=1,
        )
    fig.add_trace(
        go.Scatter(
            x=ct["time"],
            y=fit,
            marker=dict(color="gray"),
            line=dict(dash="dash"),
            mode="lines",
            name="Model",
            showlegend=True,
        ),
        row=1,
        col=1,
    )
    t1, t2 = 0, 5
    voa = fit_max_growth_rate(oa, t1, t2, plot=False)
    oa = oa[(oa["time"] >= t1) & (oa["time"] <= t2)]
    n0 = get_n0(oa)
    fit = [n0 * np.exp(voa * t) for t in oa["time"]]
    for c in oa.columns[1:]:
        fig.add_trace(
            go.Scatter(
                x=oa["time"],
                y=oa[c],
                marker=dict(color=colors["oa"]),
                mode="lines",
                showlegend=False,
            ),
            row=1,
            col=2,
        )
    fig.add_trace(
        go.Scatter(
            x=oa["time"],
            y=fit,
            marker=dict(color="gray"),
            line=dict(dash="dash"),
            mode="lines",
            name="Model",
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig = style_plot(fig)
    return vct, voa, fig


def fig3():
    ct, oa = get_ct_oa()
    vct, voa = fig2()[:2]
    fig = make_subplots(
        rows=1,
        cols=2,
        vertical_spacing=0.2,
        horizontal_spacing=0.12,
        row_heights=[0.3],
        subplot_titles=["Ct", "Oa"],
    )
    Kmct = fit_Km(ct, vct, 15)
    n0 = get_n0(ct)
    q = get_yield(ct, 15)
    y = odeint(
        monod,
        [n0, 15],
        ct["time"].to_numpy(),
        args=(vct, Kmct, q),
    )

    for c in ct.columns[1:]:
        fig.add_trace(
            go.Scatter(
                x=ct["time"],
                y=ct[c],
                marker=dict(color=colors["ct"]),
                mode="lines",
                showlegend=False,
            ),
            row=1,
            col=1,
        )
    fig.add_trace(
        go.Scatter(
            x=ct["time"],
            y=y[:, 0],
            marker=dict(color="gray"),
            line=dict(dash="dash"),
            mode="lines",
            name="Model",
            showlegend=True,
        ),
        row=1,
        col=1,
    )
    Kmoa = fit_Km(oa, voa, 15)
    n0 = get_n0(oa)
    q = get_yield(oa, 15)
    y = odeint(
        monod,
        [n0, 15],
        oa["time"].to_numpy(),
        args=(voa, Kmoa, q),
    )
    for c in oa.columns[1:]:
        fig.add_trace(
            go.Scatter(
                x=oa["time"],
                y=oa[c],
                marker=dict(color=colors["oa"]),
                mode="lines",
                showlegend=False,
            ),
            row=1,
            col=2,
        )
    fig.add_trace(
        go.Scatter(
            x=oa["time"],
            y=y[:, 0],
            marker=dict(color="gray"),
            line=dict(dash="dash"),
            mode="lines",
            name="Model",
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig = style_plot(fig)
    return Kmct, Kmoa, fig
