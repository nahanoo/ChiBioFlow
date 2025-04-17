import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from glob import glob
import pandas as pd
from fitting import *
from experiment import *

colors = {"ct": "#7570B3", "oa": "#D95F02", "total": "#565e91"}

w = 400
h = 250

r = np.array([0.37, 0.47])
K = np.array([5.1, 0.82])
q = np.array([0.024, 0.018])
M = 15
D = 0.07
N0 = [0.07, 0.07]
labels = ["Ct", "Oa"]


def model(y, t):
    R = y[0]
    N1 = y[1]
    N2 = y[2]
    dR = (
        D * M
        - D * R
        - N1 / q[0] * r[0] * R / (R + K[0])
        - N2 / q[1] * r[1] * R / (R + K[1])
    )
    dN1 = r[0] * R / (K[0] + R) * N1 - D * N1
    dN2 = r[1] * R / (K[1] + R) * N2 - D * N2
    return [dR, dN1, dN2]


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


def fig4():
    fig = go.Figure()
    t1 = 64
    xs = np.linspace(0, t1, t1 * 10)
    Y = odeint(model, [M, N0[0], N0[1]], xs)
    ct = Y[:, 1]
    oa = Y[:, 2]
    total = ct + oa

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ct,
            marker=dict(color=colors["ct"]),
            mode="lines",
            showlegend=True,
            name="Ct",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=oa,
            marker=dict(color=colors["oa"]),
            mode="lines",
            showlegend=True,
            name="Oa",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=total,
            marker=dict(color="gray"),
            mode="lines",
            showlegend=True,
            name="total",
        )
    )
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    fig = style_plot(fig)
    fig.show()


def fig5():
    (
        M0,
        M1,
        M2,
    ) = parse_data()[1:4]

    fig = make_subplots(
        rows=1,
        cols=2,
        vertical_spacing=0.2,
        horizontal_spacing=0.12,
        row_heights=[0.3],
        subplot_titles=["OD", "mCherry"],
    )
    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["od_measured"],
            marker=dict(color=colors["total"]),
            mode="lines",
            showlegend=True,
            name="Community 1",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["od_measured"],
            marker=dict(color=colors["total"]),
            mode="lines",
            showlegend=True,
            name="Community 2",
            line=dict(dash="dot"),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["od_measured"],
            marker=dict(color=colors["total"]),
            mode="lines",
            showlegend=True,
            name="Community 3",
            line=dict(dash="dash"),
        ),
        row=1,
        col=1,
    )
    t1 = 64
    xs = np.linspace(0, t1, t1 * 10)
    N0 = [0.24 / 2, 0.24 / 2]
    Y = odeint(model, [11.5, N0[0], N0[1]], xs)
    ct = Y[:, 1]
    oa = Y[:, 2]
    total = ct + oa
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=total,
            marker=dict(color="gray"),
            mode="lines",
            showlegend=True,
            name="Model",
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            mode="lines",
            showlegend=True,
            name="Community 1",
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            mode="lines",
            showlegend=True,
            name="Community 3",
            line=dict(dash="dot"),
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            mode="lines",
            showlegend=True,
            name="Community 3",
            line=dict(dash="dash"),
        ),
        row=1,
        col=2,
    )
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    fig = style_plot(fig, marker_size=2)
    fig.show()


def fig6():
    N0 = [0.24 / 2, 0.24 / 2]
    t1 = 64
    xs = np.linspace(0, t1, t1 * 10)
    Y = odeint(model, [11.5, N0[0], N0[1]], xs)
    ct, oa, total = Y[:, 1], Y[:, 2], Y[:, 1] + Y[:, 2]
    ct_rel, oa_rel = ct / total, oa / total
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ct_rel,
            marker=dict(color=colors["ct"]),
            mode="lines",
            showlegend=True,
            name="Ct",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=oa_rel,
            marker=dict(color=colors["oa"]),
            mode="lines",
            showlegend=True,
            name="Oa",
        )
    )
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="Relative abundance")
    fig = style_plot(fig)
    fig.show()


def fig7():
    fig = make_subplots(
        shared_yaxes=True,
        rows=1,
        cols=3,
        vertical_spacing=0.2,
        horizontal_spacing=0.12,
        row_heights=[0.3],
        subplot_titles=["Community 1", "Community 2", "Community 3"],
    )
    (ct_com1, oa_com1, ct_com2, oa_com2, ct_com3, oa_com3) = parse_data()[-6:]

    fig.add_trace(
        go.Scatter(
            x=ct_com1["sample_time"],
            y=ct_com1["average"],
            marker=dict(color=colors["ct"]),
            error_y=dict(type="data", array=ct_com1["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=ct_com2["sample_time"],
            y=ct_com2["average"],
            marker=dict(color=colors["ct"]),
            error_y=dict(type="data", array=ct_com2["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=ct_com3["sample_time"],
            y=ct_com3["average"],
            marker=dict(color=colors["ct"]),
            error_y=dict(type="data", array=ct_com3["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=1,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=oa_com1["sample_time"],
            y=oa_com1["average"],
            marker=dict(color=colors["oa"]),
            error_y=dict(type="data", array=oa_com1["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=oa_com2["sample_time"],
            y=oa_com2["average"],
            marker=dict(color=colors["oa"]),
            error_y=dict(type="data", array=oa_com2["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=oa_com3["sample_time"],
            y=oa_com3["average"],
            marker=dict(color=colors["oa"]),
            error_y=dict(type="data", array=oa_com3["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=1,
        col=3,
    )
    fig["layout"]["yaxis1"]["type"] = "log"
    fig["layout"]["yaxis2"]["type"] = "log"
    fig["layout"]["yaxis3"]["type"] = "log"

    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="CFUs/mL")
    fig = style_plot(fig, marker_size=2)
    fig.show()


def fig8():
    fig = make_subplots(
        shared_yaxes=True,
        rows=1,
        cols=3,
        vertical_spacing=0.2,
        horizontal_spacing=0.12,
        row_heights=[0.3],
        subplot_titles=["Community 1", "Community 2", "Community 3"],
    )
    (ct_com1, oa_com1, ct_com2, oa_com2, ct_com3, oa_com3) = parse_data()[-6:]

    fig.add_trace(
        go.Scatter(
            x=ct_com1["sample_time"],
            y=ct_com1["composition"],
            marker=dict(color=colors["ct"]),
            name="Ct",
            showlegend=True,
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=ct_com2["sample_time"],
            y=ct_com2["composition"],
            marker=dict(color=colors["ct"]),
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=ct_com3["sample_time"],
            y=ct_com3["composition"],
            marker=dict(color=colors["ct"]),
            showlegend=False,
        ),
        row=1,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=oa_com1["sample_time"],
            y=oa_com1["composition"],
            marker=dict(color=colors["oa"]),
            showlegend=True,
            name="Oa",
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=oa_com2["sample_time"],
            y=oa_com2["composition"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=oa_com3["sample_time"],
            y=oa_com3["composition"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=1,
        col=3,
    )

    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="CFUs/mL")
    fig = style_plot(fig, marker_size=2)
    fig.show()


def fig9():
    N0 = [0.0276, 0.0924]
    t1 = 64
    xs = np.linspace(0, t1, t1 * 10)
    Y = odeint(model, [11.5, N0[0], N0[1]], xs)
    ct, oa, total = Y[:, 1], Y[:, 2], Y[:, 1] + Y[:, 2]
    ct_rel, oa_rel = ct / total, oa / total
    fig = make_subplots(
        shared_yaxes=True,
        rows=1,
        cols=3,
        vertical_spacing=0.2,
        horizontal_spacing=0.12,
        row_heights=[0.3],
        subplot_titles=["Community 1", "Community 2", "Community 3"],
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ct_rel,
            marker=dict(color=colors["ct"]),
            mode="lines",
            showlegend=False,
            name="Ct",
        ),
        col=[1, 2, 3],
        row=[1, 1, 1],
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=oa_rel,
            marker=dict(color=colors["oa"]),
            mode="lines",
            showlegend=False,
            name="Oa",
        ),
        col=[1, 2, 3],
        row=[1, 1, 1],
    )
    (ct_com1, oa_com1, ct_com2, oa_com2, ct_com3, oa_com3) = parse_data()[-6:]

    fig.add_trace(
        go.Scatter(
            x=ct_com1["sample_time"],
            y=ct_com1["composition"],
            marker=dict(color=colors["ct"]),
            mode="markers",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=ct_com2["sample_time"],
            y=ct_com2["composition"],
            marker=dict(color=colors["ct"]),
            showlegend=False,
            mode="markers",
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=ct_com3["sample_time"],
            y=ct_com3["composition"],
            marker=dict(color=colors["ct"]),
            showlegend=True,
            mode="markers",
            name="Ct",
        ),
        row=1,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=oa_com1["sample_time"],
            y=oa_com1["composition"],
            marker=dict(color=colors["oa"]),
            showlegend=True,
            mode="markers",
            name="Oa",
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=oa_com2["sample_time"],
            y=oa_com2["composition"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
            mode="markers",
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=oa_com3["sample_time"],
            y=oa_com3["composition"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
            mode="markers",
        ),
        row=1,
        col=3,
    )
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="Relative abundance")
    fig = style_plot(fig)
    fig.show()


def fig10():
    global D
    D = 0.11
    N0 = [0.24 / 2, 0.24 / 2]
    t1 = 64
    xs = np.linspace(0, t1, t1 * 10)
    Y = odeint(model, [11.5, N0[0], N0[1]], xs)
    ct, oa, total = Y[:, 1], Y[:, 2], Y[:, 1] + Y[:, 2]
    ct_rel, oa_rel = ct / total, oa / total
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ct_rel,
            marker=dict(color=colors["ct"]),
            mode="lines",
            showlegend=True,
            name="Ct",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=oa_rel,
            marker=dict(color=colors["oa"]),
            mode="lines",
            showlegend=True,
            name="Oa",
        )
    )
    fig["layout"]["width"], fig["layout"]["height"] = w * 2.5, h * 1.5
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="Relative abundance")
    fig = style_plot(fig)
    fig.show()
