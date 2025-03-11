import numpy as np
from scipy.integrate import odeint
import pandas as pd
from style import *
import plotly.graph_objects as go

df = dict(pd.read_csv("parameters.csv"))
params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
p = params

width, height = 300, 250
lm = 10
bm = 10
tm = 10
rm = 10
font_size = 8
line_thickness = 1.2
xs = np.linspace(0, 1000, 2000)


def ct_mono(y, t):
    Ct, R = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    dCt = JCt * Ct - p["D"] * Ct
    dR = -JCt * Ct / p["q1_1"] - p["D"] * R + p["D"] * p["M1"]
    return dCt, dR


def oa_mono(y, t):
    Oa, R = y
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
    return dOa, dR


def plot_ct_mono():
    Y = odeint(ct_mono, [p["N01"], p["M1"]], xs)
    Ct, R = Y[-1]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 0],
            marker=dict(color=colors["ct"]),
            opacity=0.5,
            name="<i>Ct</i> monoculture",
        )
    )
    return fig


def plot_oa_mono():
    Y = odeint(oa_mono, [p["N02"], p["M1"]], xs)
    Oa, R = Y[-1]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 0],
            marker=dict(color=colors["oa"]),
            opacity=0.5,
            name="<i>Oa</i> monoculture",
        )
    )
    return fig


def competition(y, t):
    Ct, Oa, R = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    return dCt, dOa, dR


def plot_competition():
    Y = odeint(competition, [p["N01"], p["N02"], p["M1"]], xs)
    Ct, Oa, R = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    ct_mono_trace = plot_ct_mono()
    oa_mono_trace = plot_oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="Abundance [OD]")
    fig.update_layout(width=width, height=height)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/competition.svg")


def thiamine_supply(y, t):
    Ct, Oa, R, T = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    dT = -JOa * Oa / p["q2_3"] - p["D"] * T + p["M3"] * p["D"]
    return dCt, dOa, dR, dT


def plot_thiamine_supply():
    Y = odeint(thiamine_supply, [p["N01"], p["N02"], p["M1"], p["M3"]], xs)
    Ct, Oa, R, T = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    ct_mono_trace = plot_ct_mono()
    oa_mono_trace = plot_oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="Abundance [OD]")
    fig.update_layout(width=width, height=height)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/thiamine_supply.svg")


def mutual_cf(y, t):
    Ct, Oa, R, T = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    dT = p["a1_3"] * Ct / p["q1_3"] - JOa * Oa / p["q2_3"] - p["D"] * T
    return dCt, dOa, dR, dT


def plot_mutual_cf():
    Y = odeint(mutual_cf, [p["N01"], p["N02"], p["M1"], 0], xs)
    Ct, Oa, R, T = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    ct_mono_trace = plot_ct_mono()
    oa_mono_trace = plot_oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/mutual_cf.svg")


def niche_creation(y, t):
    Ct, Oa, R, M = y
    JCtR = p["v1_1"] * R / (R + p["K1_1"])
    JCtM = p["v1_2"] * M / (M + p["K1_2"])
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dCt = JCtR * Ct + JCtM * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCtR * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
    dM = p["a2_2"] * Oa / p["q2_2"] - JCtM * Ct / p["q1_2"] - p["D"] * M
    return dCt, dOa, dR, dM


def plot_niche_creation():
    Y = odeint(niche_creation, [p["N01"], p["N02"], p["M1"], 0], xs)
    Ct, Oa, R, M = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "M", M)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    Ct, Oa, R, M = Y[-1]
    ct_mono_trace = plot_ct_mono()
    oa_mono_trace = plot_oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/niche_creation.svg")


def niche_creation_cf(y, t):
    Ct, Oa, R, T, M = y
    JCtR = p["v1_1"] * R / (R + p["K1_1"])
    JCtM = p["v1_2"] * M / (M + p["K1_2"])
    JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
    dCt = JCtR * Ct + JCtM * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCtR * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
    dT = p["a1_3"] * Ct / p["q1_3"] - JOa * Oa / p["q2_3"] - p["D"] * T
    dM = p["a2_2"] * Oa / p["q2_2"] - JCtM * Ct / p["q1_2"] - p["D"] * M
    return dCt, dOa, dR, dT, dM


def plot_niche_creation_cf():
    Y = odeint(niche_creation_cf, [p["N01"], p["N02"], p["M1"], 0, 0], xs)
    Ct, Oa, R, T, M = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T, "M", M)
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    ct_mono_trace = plot_ct_mono()
    oa_mono_trace = plot_oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    Ct, Oa, R, T, M = Y[-1]
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/niche_creation_cf.svg")
