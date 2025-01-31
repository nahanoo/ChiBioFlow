import numpy as np
from scipy.integrate import odeint
import pandas as pd
from style import *
import plotly.graph_objects as go

df = dict(pd.read_csv("parameters.csv"))
params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
p = params

width, height = 2 * 300, 2 * 250
xs = np.linspace(0, 2000, 2000)


def ct_mono():
    def model(y, t):
        Ct, R = y
        JCt = p["v1_1"] * R / (R + p["K1_1"])
        dCt = JCt * Ct - p["D"] * Ct
        dR = -JCt * Ct / p["q1_1"] - p["D"] * R + p["D"] * p["M1"]
        return dCt, dR

    Y = odeint(model, [p["N01"], p["M1"]], xs)
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


def oa_mono():
    def model(y, t):
        Oa, R = y
        JOa = p["v2_1"] * R / (R + p["K2_1"])
        dOa = JOa * Oa - p["D"] * Oa
        dR = -JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
        return dOa, dR

    Y = odeint(model, [p["N02"], p["M1"]], xs)
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


def competition():
    def model(y, t):
        Ct, Oa, R = y
        JCt = p["v1_1"] * R / (R + p["K1_1"])
        JOa = p["v2_1"] * R / (R + p["K2_1"])
        dCt = JCt * Ct - p["D"] * Ct
        dOa = JOa * Oa - p["D"] * Oa
        dR = (
            -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
        )
        return dCt, dOa, dR

    Y = odeint(model, [p["N01"], p["N02"], p["M1"]], xs)
    Ct, Oa, R = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    ct_mono_trace = ct_mono()
    oa_mono_trace = oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height)
    fig = style_plot(fig)
    fig.write_image("plots/simulations/competition.svg")


def model(y, t):
    Ct, Oa, R, T = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    dT = -JOa * Oa / p["q2_3"] - p["D"] * T + p["M3"] * p["D"] + JCt * Ct / p["a1_3"]
    return dCt, dOa, dR, dT


Ms = np.linspace(0, 7.5, 100)
Ts = np.linspace(0, 10, 100)
ratios = []
for M in Ms:
    alpha_ratios = []
    for T in Ts:
        p["M1"] = M
        p["M3"] = T
        Y = odeint(model, [p["N01"], p["N02"], p["M1"], p["M3"]], xs)
        Ct, Oa, R, T = Y[-1]
        if (Ct < 1e-3) & (Oa < 1e-3):
            alpha_ratios.append(None)
        elif (Ct < 1e-3) & (Oa > 1e-3):
            alpha_ratios.append(1)
        elif (Ct > 1e-3) & (Oa < 1e-3):
            alpha_ratios.append(0)
        else:
            alpha_ratios.append(Oa / (Ct + Oa))
    ratios.append(alpha_ratios)
custom_colorscale = [
    [0, "blue"],  # Min value -> Blue
    [0.5, "white"],  # Mid value (zero) -> White
    [1, "red"],  # Max value -> Red
]
fig = go.Figure()
fig.add_trace(
    go.Contour(
        z=np.array(ratios),
        x=Ts,
        y=Ms,
        colorscale=custom_colorscale,
        ncontours=50,
        contours=dict(showlines=False),
        colorbar=dict(
            title=dict(
                text="<i>Oa</i> fraction",
                side="right",
            ),
            len=0.5,
        ),
    )
)
fig.write_image("plots/simulations/supp_thiamine_ratios.svg")


def mutual_cf():
    def model(y, t):
        Ct, Oa, R, T = y
        JCt = p["v1_1"] * R / (R + p["K1_1"])
        JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
        dCt = JCt * Ct - p["D"] * Ct
        dOa = JOa * Oa - p["D"] * Oa
        dR = (
            -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
        )
        dT = JCt * Ct / p["a1_3"] - JOa * Oa * p["q2_3"] - p["D"] * T
        return dCt, dOa, dR, dT

    Y = odeint(model, [p["N01"], p["N02"], p["M1"], 0], xs)
    Ct, Oa, R, T = Y[-1]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    ct_mono_trace = ct_mono()
    oa_mono_trace = oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height)
    fig = style_plot(fig)
    fig.write_image("plots/simulations/mutual_cf.svg")


def cf_mutual_species_ratio():
    def model(y, t):
        Ct, Oa, R, T = y
        JCt = p["v1_1"] * R / (R + p["K1_1"])
        JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
        dCt = JCt * Ct - p["D"] * Ct
        dOa = JOa * Oa - p["D"] * Oa
        dR = (
            -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
        )
        dT = JCt * Ct * p["q1_3"] - JOa * Oa * p["q2_3"] - p["D"] * T
        return dCt, dOa, dR, dT

    Ms = np.linspace(0, 7.5, 100)
    alphas = np.linspace(0, 5, 100)
    ratios = []
    for M in Ms:
        alpha_ratios = []
        for alpha in alphas:
            p["M1"] = M
            p["q1_3"] = alpha
            Y = odeint(model, [p["N01"], p["N02"], p["M1"], 0], xs)
            Ct, Oa, R, T = Y[-1]
            if (Ct < 1e-3) & (Oa < 1e-3):
                alpha_ratios.append(None)
            else:
                alpha_ratios.append(Oa / Ct)
        ratios.append(alpha_ratios)
    custom_colorscale = [
        [0, "blue"],  # Min value -> Blue
        [0.5, "white"],  # Mid value (zero) -> White
        [1, "red"],  # Max value -> Red
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=np.array(ratios),
            x=alphas,
            y=Ms,
            colorscale=custom_colorscale,
            ncontours=50,
            contours=dict(showlines=False),
            colorbar=dict(
                title=dict(
                    text="<i>Oa</i> fraction",
                    side="right",
                ),
                len=0.5,
            ),
        )
    )
    fig.write_image("plots/simulations/mutual_cf_ratios.svg")


def niche_creation():
    def model(y, t):
        Ct, Oa, R, M = y
        JCtR = p["v1_1"] * R / (R + p["K1_1"])
        JCtM = p["v1_2"] * M / (M + p["K1_2"])
        JOa = p["v2_1"] * R / (R + p["K2_1"])
        dCt = JCtR * Ct + JCtM * Ct - p["D"] * Ct
        dOa = JOa * Oa - p["D"] * Oa
        dR = (
            -JCtR * Ct / p["q1_1"]
            - JOa * Oa / p["q2_1"]
            - p["D"] * R
            + p["D"] * p["M1"]
        )
        dM = JOa * Oa / p["a2_2"] - JCtM * Ct / p["q1_2"] - p["D"] * M
        return dCt, dOa, dR, dM

    Y = odeint(model, [p["N01"], p["N02"], p["M1"], 0], xs)
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
    ct_mono_trace = ct_mono()
    oa_mono_trace = oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height)
    fig = style_plot(fig)
    fig.write_image("plots/simulations/niche_creation.svg")


def niche_creation_cf():
    def model(y, t):
        Ct, Oa, R, T, M = y
        JCtR = p["v1_1"] * R / (R + p["K1_1"])
        JCtM = p["v1_2"] * M / (M + p["K1_2"])
        JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
        dCt = JCtR * Ct + JCtM * Ct - p["D"] * Ct
        dOa = JOa * Oa - p["D"] * Oa
        dR = (
            -JCtR * Ct / p["q1_1"]
            - JOa * Oa / p["q2_1"]
            - p["D"] * R
            + p["D"] * p["M1"]
        )
        dT = (
            JCtM * Ct / p["a1_3"]
            + JCtR * Ct / p["a1_3"]
            - JOa * Oa / p["q2_3"]
            - p["D"] * T
        )
        dM = JOa * Oa / p["a2_2"] - JCtM * Ct / p["q1_2"] - p["D"] * M
        return dCt, dOa, dR, dT, dM

    Y = odeint(model, [p["N01"], p["N02"], p["M1"], 0, 0], xs)
    Ct, Oa, R, T, M = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T, "M", M)
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    ct_mono_trace = ct_mono()
    oa_mono_trace = oa_mono()
    fig.add_trace(ct_mono_trace.data[0])
    fig.add_trace(oa_mono_trace.data[0])
    Ct, Oa, R, T, M = Y[-1]
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height)
    fig = style_plot(fig)
    fig.write_image("plots/simulations/niche_creation_cf.svg")
    return R, M
