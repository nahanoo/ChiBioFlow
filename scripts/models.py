import numpy as np
from scipy.integrate import odeint
import pandas as pd
from style import *
import plotly.graph_objects as go
import plotly.io as pio


def parse_params():
    df = dict(pd.read_csv("parameters.csv"))
    params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
    p = params
    return p


# width, height = 300, 250
lm = 10
bm = 10
tm = 10
rm = 10
font_size = 8
line_thickness = 1.2
xs = np.linspace(0, 5000, 5000 * 6)


def ct_mono(y, t, p):
    Ct, R = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    dCt = JCt * Ct - p["D"] * Ct
    dR = -JCt * Ct / p["q1_1"] - p["D"] * R + p["D"] * p["M1"]
    return dCt, dR


def oa_mono(y, t, p):
    Oa, R = y
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
    return dOa, dR


def plot_ct_mono():
    p = parse_params()
    Y = odeint(ct_mono, [p["N01"], p["M1"]], xs, args=(p,))
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
    p = parse_params()
    Y = odeint(oa_mono, [p["N02"], p["M1"]], xs, args=(p,))
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


def competition(y, t, p):
    Ct, Oa, R = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    return dCt, dOa, dR


def plot_competition():
    p = parse_params()
    p["D"] = 0.15
    xs = np.linspace(0, 250, 1000)
    Y = odeint(competition, [p["N01"], p["N02"], p["M1"]], xs, args=(p,))
    Ct, Oa, R = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="Abundance [OD]")
    fig.update_layout(
        width=width,
        height=height,
        # xaxis=dict(range=[0, 72], dtick=24),
        yaxis=dict(range=[0, 0.3], dtick=0.1),
        showlegend=False,
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/dynamics/fig3b.svg")


def thiamine_supply(y, t, p):
    Ct, Oa, R, T = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    dT = -JOa * Oa / p["q2_3"] - p["D"] * T + p["M3"] * p["D"]
    return dCt, dOa, dR, dT


def plot_thiamine_supply():
    p = parse_params()
    Y = odeint(thiamine_supply, [p["N01"], p["N02"], p["M1"], p["M3"]], xs, args=(p,))
    Ct, Oa, R, T = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
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
    fig.write_image("plots/simulations/dynamics/thiamine_supply.svg")


def plot_thiamine_supply_oa_batch():
    xs = np.linspace(0, 72, 2000)
    dash = ["dash", "dot", "dashdot"]
    names = ["0 nM thiamine", "1 nM thiamine", "10 nM thiamine"]
    ts = [0, 1, 10]
    p = parse_params()
    p["v1_1"] = 0
    p["D"] = 0
    fig = go.Figure()
    for i, t in enumerate(ts):
        p["M3"] = t
        Y = odeint(
            thiamine_supply, [p["N01"], p["N02"], p["M1"], p["M3"]], xs, args=(p,)
        )
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=Y[:, 1],
                name=names[i],
                marker=dict(color=colors["oa"]),
                line=dict(dash=dash[i]),
            )
        )
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
    fig.write_image("plots/simulations/parameterization/thiamine_supply_oa_batch.pdf")


def mutual_cf(y, t, p):
    Ct, Oa, R, T = y
    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])
    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    dT = JCt * Ct / p["q1_3"] - JOa * Oa / p["q2_3"] - p["D"] * T
    return dCt, dOa, dR, dT


def plot_mutual_cf():
    xs = np.linspace(0, 250, 2000)
    p = parse_params()
    p["N02"] = 0
    Y = odeint(mutual_cf, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
    Ct, Oa, R, T = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
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
    fig.write_image("plots/simulations/dynamics/mutual_cf.svg")


def plot_thiamine_production():
    xs = np.linspace(0, 200, 2000)
    p = parse_params()
    p["v2_1"] = 0
    Y = odeint(mutual_cf, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
    Ct, Oa, R, T = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 3],
            name="Thiamine<br>concentration",
            marker=dict(color=colors["blue"]),
        )
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="Concentration [nM]")
    fig.update_layout(width=width, height=height, title="Thiamine concentration")
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/parameterization/thiamine_production.pdf")


def niche_creation(y, t, p):
    Ct, Oa, R, M = y
    JCtR = p["v1_1"] * R / (R + p["K1_1"])
    JCtM = p["v1_2"] * M / (M + p["K1_2"])
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dCt = JCtR * Ct + JCtM * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCtR * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
    dM = p["a2_2"] * Oa * JOa / p["q2_2"] - JCtM * Ct / p["q1_2"] - p["D"] * M
    return dCt, dOa, dR, dM


def plot_niche_creation():
    p = parse_params()
    p["a2_2"] = 0.027
    xs = np.linspace(0, 250, 2000)
    Y = odeint(niche_creation, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
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
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="OD")
    fig.update_layout(width=width, height=height, showlegend=False)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/dynamics/fig3b.svg")


plot_niche_creation()


def plot_niche_production():
    xs = np.linspace(0, 2000, 2000)
    p = parse_params()
    p["v1_1"] = 0
    p["v1_2"] = 0
    Y = odeint(niche_creation, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
    Ct, Oa, R, M = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "M", M)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 3],
            name="Metabolite<br>concentration",
            marker=dict(color=colors["blue"]),
        )
    )
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(title="Concentration [mM]")
    fig.update_layout(width=width, height=height, title="Metabolite concentration")
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/parameterization/niche_concentration.pdf")


def niche_supply(y, t, p):
    Ct, Oa, R, M = y
    JCtR = p["v1_1"] * R / (R + p["K1_1"])
    JCtM = p["v1_2"] * M / (M + p["K1_2"])
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dCt = JCtR * Ct + JCtM * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCtR * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
    dM = -JCtM * Ct / p["q1_2"] - p["D"] * M + p["D"] * p["M2"]
    return dCt, dOa, dR, dM


def plot_niche_supply():
    p = parse_params()
    Y = odeint(niche_supply, [p["N01"], p["N02"], p["M1"], p["M2"]], xs, args=(p,))
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
    fig.write_image("plots/simulations/dynamics/niche_creation_supply.pdf")


def niche_creation_cf(y, t, p):
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
    p = parse_params()
    Y = odeint(niche_creation_cf, [p["N01"], p["N02"], p["M1"], 0, 0], xs, args=(p,))
    Ct, Oa, R, T, M = Y[-1]
    print("Ct", Ct, "Oa", Oa, "R", R, "T", T, "M", M)
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 0], name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )
    fig.add_trace(
        go.Scatter(x=xs, y=Y[:, 1], name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
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
    fig.write_image("plots/simulations/dynamics/niche_creation_cf.pdf")


def niche_creation_batch(y, t, p):
    Ct, Oa, R, M = y
    JCtR = p["v1_1"] * R / (R + p["K1_1"])
    JCtM = p["v1_2"] * M / (M + p["K1_2"])
    JOa = p["v2_1"] * R / (R + p["K2_1"])
    dCt = JCtR * Ct + JCtM * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCtR * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] - p["D"] * R + p["D"] * p["M1"]
    dM = p["a2_2"] * Oa / p["q2_2"] * JOa - JCtM * Ct / p["q1_2"] - p["D"] * M
    return dCt, dOa, dR, dM


def plot_niche_creation_batch():
    p = parse_params()
    p["a2_2"] = 0.027
    p["D"] = 0.1
    xs = np.linspace(0, 250, 1000)
    Y = odeint(niche_creation_batch, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
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
    fig.update_xaxes(title="Time [h]"), fig.update_yaxes(
        title="OD", range=[0, 0.3], dtick=0.05
    )
    fig.update_layout(width=width, height=height, showlegend=False)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/dynamics/niche_creation_batch.svg")


def caller():
    plot_competition()
    plot_thiamine_supply()
    plot_thiamine_supply_oa_batch()
    plot_mutual_cf()
    plot_thiamine_production()
    plot_niche_creation()
    plot_niche_production()
    plot_niche_supply()
    plot_niche_creation_cf()
