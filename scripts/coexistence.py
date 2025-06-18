import numpy as np
from scipy.integrate import odeint
import pandas as pd
from style import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from models import (
    competition as cp,
    thiamine_supply as ts,
    mutual_cf as mc,
    niche_creation as nc,
    niche_creation_cf as nccf,
    niche_supply as ns,
    niche_creation_batch as ncb,
)
from joblib import Parallel, delayed
import plotly.io as pio

pio.kaleido.scope.mathjax = None


def parse_params():
    df = dict(pd.read_csv("parameters.csv"))
    params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
    p = params
    return p


lm = 10
bm = 10
tm = 10
rm = 10
font_size = 8
line_thickness = 1.2
xs = np.linspace(0, 5000, 5000 * 6)

custom_colorscale = [
    [0, "#c0bedc"],
    [0.5, "white"],
    [1, "#ecaf80"],
]


def competition():
    p = parse_params()
    Ds = np.linspace(0, 0.4, 100)
    Cts = []
    Oas = []
    for D in Ds:
        p["D"] = D
        Y = odeint(cp, [p["N01"], p["N02"], p["M1"]], xs, args=(p,))
        Ct, Oa, R = Y[-1]
        Cts.append(Ct), Oas.append(Oa)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=Ds, y=Cts, name="<i>Ct</i>", marker=dict(color=colors["ct"]))
    )

    fig.add_trace(
        go.Scatter(x=Ds, y=Oas, name="<i>Oa</i>", marker=dict(color=colors["oa"]))
    )
    fig.update_xaxes(title="Dilution rate [1/h]"), fig.update_yaxes(
        title="Abundance [OD]"
    )
    fig.update_layout(
        width=width,
        height=height,
    )
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/competition.pdf")


def thiamine_supply():
    p = parse_params()
    Ts = np.linspace(0, 100, 10000)
    fig = go.Figure()
    js = p["v2_1"] * 0.1 / (0.1 + p["K2_1"]) * Ts / (Ts + p["K2_3"])
    js[-1] = 0.25
    Ts[-1] = 100
    fig.add_trace(go.Scatter(x=list(Ts), y=list(js), mode="lines"))
    fig.update_xaxes(
        title="R<sub>0,T</sub> [nM]",
        zeroline=True,
        showgrid=False,
        range=[0, 110],
        dtick=50,
    ), fig.update_yaxes(title="D [1/h]", zeroline=True, showgrid=False)
    fig.update_layout(height=150, width=120)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=40,
        buttom_margin=25,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("plots/simulations/coexistence/thiamine_supply.svg")


def mutual_cf():
    p = parse_params()
    Ds = np.linspace(0, 0.2, 50)
    alphas = 1 / np.linspace(0.002, 0.02, 50)
    zs = np.zeros((len(Ds), len(alphas)))
    for i, D in enumerate(Ds):
        p["D"] = D
        for j, alpha in enumerate(alphas):
            p["q1_3"] = alpha
            Y = odeint(mc, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
            Ct, Oa, R, T = Y[-1]
            if Ct <= 1e-6:
                Ct = 0
            if Oa <= 1e-6:
                Oa = 0
            if (Ct == 0) and (Oa == 0):
                ratio = None
            else:
                ratio = Oa / (Ct + Oa)
            zs[i, j] = ratio

    fig = go.Figure()

    fig.add_trace(
        go.Contour(
            z=zs,
            x=alphas,
            y=Ds,
            colorscale=custom_colorscale,
            ncontours=50,
            zmid=0.5,
            zmin=0,
            zmax=1,
            contours=dict(
                showlines=False,
            ),
            colorbar=dict(
                title=dict(text="<i>Oa</i> fraction", side="right", font=dict(size=8)),
                len=0.8,
                # y=0.25,
                thickness=10,
            ),
        )
    )

    fig.update_xaxes(
        title="1 / Q<sub>Ct,R</sub> nM/OD",
        zeroline=False,
    )
    fig.update_yaxes(title="Dilution rate [1/h]", zeroline=False, showgrid=False)
    fig.update_layout(height=150, width=170)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=40,
        buttom_margin=25,
        top_margin=0,
        right_margin=0,
    )
    fig.write_image("plots/simulations/coexistence/mutual_cf.svg")


def niche_creation():
    p = parse_params()
    Ds = np.linspace(0, 0.3, 100)
    alphas = np.linspace(0, 1, 100)
    zs = np.zeros((len(Ds), len(alphas)))
    for i, D in enumerate(Ds):
        p["D"] = D
        for j, alpha in enumerate(alphas):
            p["a2_2"] = alpha
            Y = odeint(nc, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
            Ct, Oa, R, M = Y[-1]
            ratio = Oa / (Ct + Oa)
            if ratio <= 1e-5:
                ratio = 0
            zs[i, j] = ratio

    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=zs,
            x=alphas,
            y=Ds,
            colorscale=custom_colorscale,
            ncontours=50,
            zmid=0.5,
            zmin=0,
            zmax=1,
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
    fig.update_xaxes(title="Metabolite leakage rate [1/h]"), fig.update_yaxes(
        title="Dilution rate [1/h]"
    )
    fig.update_layout(
        height=height,
        width=width,
        xaxis=dict(showgrid=False, ticks="outside"),
        yaxis=dict(showgrid=False, ticks="outside"),
    )
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/niche_creation.pdf")


def niche_creation_cf():
    p = parse_params()
    alphas = np.linspace(0, 1, 10)
    zs = np.zeros((len(alphas), len(alphas)))
    for i, alpha in enumerate(alphas):
        p["a1_3"] = alpha
        for j, beta in enumerate(alphas):
            p["a2_2"] = beta
            Y = odeint(nccf, [p["N01"], p["N02"], p["M1"], 0, 0], xs, args=(p,))
            Ct, Oa, R, T, M = Y[-1]
            ratio = Oa / (Ct + Oa)
            if ratio <= 1e-5:
                ratio = 0
            zs[i, j] = ratio

    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=zs,
            x=alphas,
            y=alphas,
            colorscale=custom_colorscale,
            zmid=0.5,
            zmin=0,
            zmax=1,
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
    fig.update_xaxes(title="Metabolite leakage rate [1/h]"), fig.update_yaxes(
        title="Thiamine leakage rate [1/h]"
    )
    fig.update_layout(
        height=height,
        width=width,
        xaxis=dict(showgrid=False, ticks="outside"),
        yaxis=dict(showgrid=False, ticks="outside"),
    )
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/niche_creation_cf.pdf")


def misleading_interaction():
    fig = make_subplots(
        rows=1,
        cols=4,
        horizontal_spacing=0.1,
    )
    figj = make_subplots(
        rows=1,
        cols=3,
        horizontal_spacing=0.1,
    )
    p = parse_params()
    a = 0.027
    p["D"] = 0
    p["N02"] = 0
    p["a2_2"] = a
    xs = np.linspace(0, 24, 1000)
    Y = odeint(nc, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
    R = Y[:, 2]
    JCt = p["v1_1"] * R / (p["K1_1"] + R)
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 0],
            name="Ct",
            line=dict(color=colors["ct"], shape="spline"),
        ),
        row=1,
        col=1,
    )
    figj.add_trace(
        go.Scatter(
            x=xs,
            y=JCt,
            name="Ct",
            line=dict(color=colors["ct"], shape="spline"),
        ),
        row=1,
        col=1,
    )

    p = parse_params()
    p["D"] = 0
    p["N01"] = 0
    p["a2_2"] = a

    Y = odeint(nc, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
    R = Y[:, 2]
    JOa = p["v2_1"] * R / (p["K2_1"] + R)
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 1],
            name="Oa",
            line=dict(color=colors["oa"], shape="spline"),
        ),
        row=1,
        col=2,
    )
    figj.add_trace(
        go.Scatter(
            x=xs,
            y=JOa,
            name="Oa",
            line=dict(color=colors["oa"], shape="spline"),
        ),
        row=1,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 3],
            name="M",
            line=dict(color="black", shape="spline"),
        ),
        row=1,
        col=4,
    )
    p = parse_params()
    p["D"] = 0
    p["a2_2"] = a

    Y = odeint(nc, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
    R = Y[:, 2]
    JCt = p["v1_1"] * R / (p["K1_1"] + R)
    JOa = p["v2_1"] * R / (p["K2_1"] + R)

    M = Y[:, 3]
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 0],
            name="Ct",
            mode="lines",
            line=dict(color=colors["ct"], shape="spline"),
            showlegend=False,
        ),
        row=1,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 1],
            name="Oa",
            line=dict(color=colors["oa"], shape="spline"),
            showlegend=False,
        ),
        row=1,
        col=3,
    )

    figj.add_trace(
        go.Scatter(
            x=xs,
            y=JOa,
            name="Oa",
            line=dict(color=colors["oa"], shape="spline"),
            showlegend=False,
        ),
        row=1,
        col=3,
    )
    figj.add_trace(
        go.Scatter(
            x=xs,
            y=JCt,
            name="Ct",
            line=dict(color=colors["ct"], shape="spline"),
            showlegend=False,
        ),
        row=1,
        col=3,
    )

    fig.for_each_xaxis(lambda x: x.update(range=[0, 24], dtick=12))
    fig.for_each_yaxis(lambda y: y.update(range=[0, 0.4], dtick=0.2))
    fig.update_layout(width=width * 2, height=height, showlegend=False)
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=lm,
        buttom_margin=0,
        top_margin=10,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/batch_failure.svg")

    figj.for_each_xaxis(lambda x: x.update(range=[0, 24], dtick=12))
    figj.for_each_yaxis(lambda y: y.update(range=[0, 0.4], dtick=0.2))
    figj.update_layout(width=width * 2, height=height, showlegend=False)
    figj = style_plot(
        figj,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=lm,
        buttom_margin=0,
        top_margin=10,
        right_margin=rm,
    )
    figj.write_image("plots/simulations/coexistence/batch_failure_j.svg")


def differences():
    Ds = np.linspace(0, 0.2, 100)
    p = parse_params()
    JCts = []
    JOas = []
    JCts_diff = []
    for D in Ds:
        p["D"] = D
        Y = odeint(cp, [p["N01"], p["N02"], p["M1"]], xs, args=(p,))
        R = Y[:, 2][-1]
        JCt = p["v1_1"] * R / (p["K1_1"] + R)
        JOa = p["v2_1"] * R / (p["K2_1"] + R)
        JCts.append(JCt)
        JOas.append(JOa)
        JCts_diff.append(D - JCt)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=Ds, y=JCts, name="Ct", marker=dict(color=colors["ct"]), mode="lines"
        )
    )
    fig.add_trace(
        go.Scatter(
            x=Ds, y=JOas, name="Oa", marker=dict(color=colors["oa"]), mode="lines"
        )
    )
    fig.add_trace(
        go.Scatter(
            x=Ds,
            y=JCts_diff,
            name="D - J<sub>Ct</sub>",
            line=dict(dash="dot"),
            marker=dict(color=colors["ct"]),
            mode="lines",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Dilution rate [1/h]"),
        yaxis=dict(title="J [1/h]"),
        showlegend=False,
        width=width,
        height=height,
    )
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=10,
        buttom_margin=10,
        top_margin=10,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/differences_J.svg")


def realized_J():
    Kms = np.linspace(1e-3, 1, 1000)
    Rs = np.linspace(1e-3, 1, 1000)
    zs = np.zeros((len(Rs), len(Kms)))
    for i, R in enumerate(Rs):
        for j, Km in enumerate(Kms):
            zs[i, j] = 0.2 * R / (R + Km)
    D = 0.1
    p = parse_params()
    R_star = -D * p["K2_1"] / (D - p["v2_1"])
    J_diff = D - p["v1_1"] * R_star / (p["K1_1"] + R_star)

    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=zs,
            x=Kms,
            y=Rs,
            colorscale=custom_colorscale,
            ncontours=20,
            zmin=0,
            zmax=0.2,
            zmid=0.1,
            contours=dict(showlines=False),
            colorbar=dict(
                title=dict(text="J", side="right", font=dict(size=8)),
                len=0.6,
                y=0.25,
                thickness=10,
            ),
        )
    )

    fig.add_trace(
        go.Contour(
            z=zs,
            x=Kms,
            y=Rs,
            showscale=False,
            contours=dict(start=0.07, end=0.07, size=0.1, coloring="none"),
            line=dict(color="black"),
            showlegend=False,
        )
    )

    fig.update_layout(
        xaxis=dict(title="K<sub>M</sub> [mM]", zeroline=False, type="log", dtick="1"),
        yaxis=dict(title="M [mM]", zeroline=False, type="log", dtick="1"),
        height=150,
        width=170,
    )

    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=lm,
        buttom_margin=25,
        top_margin=5,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/realized_J.svg")


def km_dataset():
    color_dict = {
        "arabinose": "#1f77b4",  # blue
        "fructose": "#ff7f0e",  # orange
        "glucose": "#2ca02c",  # green
        "lactate": "#d62728",  # red
        "lactose": "#9467bd",  # purple
        "maltose": "#8c564b",  # brown
    }
    symbol_dict = {
        "Achromobacter sp.": "pentagon",
        "Escherichia coli": "square",
        "Marine coryneform baterium": "diamond",
        "Pseudomonas sp.": "cross",
        "Spirillum sp.": "x",
        "Streptococcus mutans": "triangle-up",
        "Streptococcus sanguis": "triangle-down",
        "Vibrio sp.": "star",
    }
    df = pd.read_excel("km_dataset.xlsx")
    df = df[df["Eukaryote"] != True]

    fig = go.Figure()

    # === Legend for nutrients (colors) ===
    for nutrient, color in color_dict.items():
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(color=color, symbol="circle", size=10),
                name=nutrient.capitalize(),
                legendgroup="color_legend",
                showlegend=True,
            )
        )

    # === Legend for species (symbols) ===
    for species, symbol in symbol_dict.items():
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(color="black", symbol=symbol, size=10),
                name=species,
                legendgroup="symbol_legend",
                showlegend=True,
            )
        )

    # === Actual data ===
    for _, row in df.iterrows():
        fig.add_trace(
            go.Scatter(
                x=[row["gmax (per hour)"]],
                y=[row["K common unit (uM)"]],
                mode="markers",
                marker=dict(
                    color=color_dict[row["Limiting Nutrient"]],
                    symbol=symbol_dict[row["Species"]],
                    size=10,
                ),
                showlegend=False,
            )
        )

    fig.update_layout(
        yaxis=dict(type="log", title="K [uM]"),
        xaxis=dict(title="max. growth rate [1/h]", range=[0, 2], dtick=0.5),
        showlegend=False,
        width=1 * width,
        height=1 * height,
    )
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=lm,
        buttom_margin=25,
        top_margin=5,
        right_margin=rm,
        marker_size=4,
    )
    fig.write_image("plots/experiments/km_dataset.svg")
