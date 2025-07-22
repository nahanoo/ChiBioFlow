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


def fig1d():

    p = parse_params()
    Ds = np.linspace(0, 0.2, 500)
    alphas = np.linspace(1, 100, 500)
    zs = np.zeros((len(Ds), len(alphas)))
    for i, D in enumerate(Ds):
        p["D"] = D
        for j, alpha in enumerate(alphas):
            p["M3"] = alpha
            Y = odeint(ts, [p["N01"], p["N02"], p["M1"], p["M3"]], xs, args=(p,))
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
            colorscale=colors_heatmap,
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
    fig.write_image("plots/simulations/coexistence/fig1d.svg")


def fig1c():
    p = parse_params()
    Ds = np.linspace(0, 0.3, 500)
    alphas = np.linspace(0.0002, 1, 500)
    zs = np.zeros((len(Ds), len(alphas)))
    Ts = []
    for i, D in enumerate(Ds):
        p["D"] = D
        for j, alpha in enumerate(alphas):
            p["q1_3"] = alpha
            Y = odeint(mc, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
            Ct, Oa, R, T = Y[-1]

            ratio = Oa / (Ct + Oa)
            zs[i, j] = ratio
            Ts.append(T)

    fig = go.Figure()
    from scipy.ndimage import zoom

    fig.add_trace(
        go.Contour(
            z=zs,
            x=Ts,
            y=Ds,
            colorscale=colors_heatmap,
            ncontours=50,
            zmid=0.5,
            zmin=0,
            zmax=1,
            contours=dict(showlines=False),
            colorbar=dict(
                title=dict(text="<i>Oa</i> fraction", side="right", font=dict(size=8)),
                # y=0.25,
                thickness=10,
                outlinewidth=0.5,
                outlinecolor="black",
            ),
            showscale=False,
        )
    )

    fig.update_xaxes(
        title="Thiamine concentration in chemostat [nM]",
        type="log",
        ticks="inside",
        # zeroline=False,
        # range=[0.0002, 0.02],
        # dtick=0.002,
    )
    fig.update_yaxes(
        title="Dilution rate [1/h]", zeroline=False, showgrid=False, ticks="inside"
    )
    fig.update_layout(height=150, width=150, title="Cross-feeding")
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=20,
        buttom_margin=25,
        top_margin=20,
        right_margin=10,
    )
    fig.write_image("plots/simulations/coexistence/fig1c.svg")


fig1c()


def fig3a():
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
        font_size=11,
        left_margin=lm,
        buttom_margin=0,
        top_margin=10,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/fig3a.svg")

    figj.for_each_xaxis(lambda x: x.update(range=[0, 24], dtick=12))
    figj.for_each_yaxis(lambda y: y.update(range=[0, 0.4], dtick=0.2))
    figj.update_layout(width=width * 2, height=height, showlegend=False)
    figj = style_plot(
        figj,
        font_size=11,
        left_margin=lm,
        buttom_margin=0,
        top_margin=10,
        right_margin=rm,
    )
    figj.write_image("plots/simulations/coexistence/batch_failure_j.svg")


def fig2a():
    Ds = np.linspace(0, 0.3, 1000)
    p = parse_params()
    r_stars = []
    for D in Ds:
        p["D"] = D
        Y = odeint(cp, [p["N01"], p["N02"], p["M1"]], xs, args=(p,))
        R = Y[:, 2][-1]
        r_stars.append(R)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=Ds,
            y=r_stars,
            name="r_star",
            line=dict(color="black", shape="spline"),
            mode="lines",
        )
    )
    r_equal = (p["K1_1"] * p["v2_1"] - p["K2_1"] * p["v1_1"]) / (p["v1_1"] - p["v2_1"])
    dc_oa = r_equal * p["v1_1"] / (p["K1_1"] + r_equal)
    fig.update_layout(
        xaxis=dict(title="Dilution rate [1/h]", ticks="inside", showgrid=False),
        yaxis=dict(title="Acetate concentration [mM]", ticks="inside", showgrid=False),
        showlegend=False,
        width=width,
        height=height * 1.3,
        title="Competition for acetate",
        shapes=[
            # Blue background from x=0 to x=0.1
            dict(
                type="rect",
                xref="x",
                yref="paper",
                x0=0,
                x1=dc_oa,
                y0=0,
                y1=1,
                fillcolor=colors["oa"],
                opacity=0.3,
                layer="below",
                line_width=0,
            ),
            # Red background from x=0.1 to x=0.3
            dict(
                type="rect",
                xref="x",
                yref="paper",
                x0=dc_oa,
                x1=max(Ds),
                y0=0,
                y1=1,
                fillcolor=colors["ct"],
                opacity=0.3,
                layer="below",
                line_width=0,
            ),
            dict(
                type="line",
                yref="paper",
                x0=dc_oa,
                y0=0,
                x1=dc_oa,
                y1=1,
                line=dict(color="black", width=1.5, dash="dot"),
            ),
        ],
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=20,
        buttom_margin=10,
        top_margin=20,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/fig2a.svg")


def fig2b():
    p = parse_params()
    r_equal = (p["K1_1"] * p["v2_1"] - p["K2_1"] * p["v1_1"]) / (p["v1_1"] - p["v2_1"])
    dc_oa = r_equal * p["v1_1"] / (p["K1_1"] + r_equal)
    Ds = np.linspace(0, dc_oa, 100)
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
            x=Ds,
            y=JOas,
            name="Oa",
            marker=dict(color=colors["oa"]),
            mode="lines",
            fill="tonexty",
            # fillcolor="white",
            fillcolor="rgba(117, 112, 179, 0.3)",
        )
    )
    JCts = []
    JOas = []
    JCts_diff = []
    Ds = np.linspace(dc_oa, 0.3, 100)

    for D in Ds:
        p["D"] = D
        Y = odeint(cp, [p["N01"], p["N02"], p["M1"]], xs, args=(p,))
        R = Y[:, 2][-1]
        JCt = p["v1_1"] * R / (p["K1_1"] + R)
        JOa = p["v2_1"] * R / (p["K2_1"] + R)
        JCts.append(JCt)
        JOas.append(JOa)
        JCts_diff.append(D - JCt)
    fig.add_trace(
        go.Scatter(
            x=Ds, y=JCts, name="Ct", marker=dict(color=colors["ct"]), mode="lines"
        )
    )
    fig.add_trace(
        go.Scatter(
            x=Ds,
            y=JOas,
            name="Oa",
            marker=dict(color=colors["oa"]),
            mode="lines",
            fill="tonexty",
            # fillcolor="white",
            fillcolor="rgba(217, 95, 2, 0.3)",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Dilution rate [1/h]", ticks="inside"),
        yaxis=dict(title="J [1/h]", ticks="inside"),
        showlegend=False,
        width=width,
        height=height * 1.3,
        title="Growth rate based on acetate",
        shapes=[
            dict(
                type="line",
                x0=dc_oa,
                yref="paper",
                y0=0,
                x1=dc_oa,
                y1=1,
                line=dict(color="black", width=1.5, dash="dot"),
            )
        ],
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=10,
        buttom_margin=10,
        top_margin=10,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/fig2b.svg")


def fig2c():
    p = parse_params()
    r_equal = (p["K1_1"] * p["v2_1"] - p["K2_1"] * p["v1_1"]) / (p["v1_1"] - p["v2_1"])
    dc_oa = r_equal * p["v1_1"] / (p["K1_1"] + r_equal)
    Ds = np.linspace(0, 0.3, 1000)
    p = parse_params()
    JCts = []
    JOas = []
    JCts_diff = []
    JOas_diff = []
    for D in Ds:
        p["D"] = D
        Y = odeint(cp, [p["N01"], p["N02"], p["M1"]], xs, args=(p,))
        R = Y[:, 2][-1]
        JCt = p["v1_1"] * R / (p["K1_1"] + R)
        JOa = p["v2_1"] * R / (p["K2_1"] + R)
        JCts.append(JCt)
        JOas.append(JOa)
        JCts_diff.append(D - JCt)
        JOas_diff.append(D - JOa)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=Ds,
            y=JCts_diff,
            name="Ct",
            line=dict(color=colors["ct"], shape="spline"),
            mode="lines",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=Ds,
            y=JOas_diff,
            name="Oa",
            line=dict(color=colors["oa"], shape="spline"),
            mode="lines",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Dilution rate [1/h]", ticks="inside"),
        yaxis=dict(title="J 1/h", ticks="inside"),
        showlegend=False,
        width=width,
        height=height * 1.3,
        title="Missing growth rate<br>for coexistence",
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=20,
        buttom_margin=10,
        top_margin=30,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/fig2c.svg")


def fig2e():
    Kms = np.linspace(1e-3, 1, 1000)
    Rs = np.linspace(1e-3, 1, 1000)
    zs = np.zeros((len(Rs), len(Kms)))
    for i, R in enumerate(Rs):
        for j, Km in enumerate(Kms):
            zs[i, j] = 0.4 * R / (R + Km)
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
            colorscale=[
                [0.0, "white"],  # deep blue
                [1.0, "black"],  # deep red
            ],
            ncontours=100,
            # zmin=0,
            # zmax=0.4,
            # zmid=0.2,
            contours=dict(showlines=False),
            colorbar=dict(
                title=dict(text="J", side="right", font=dict(size=11)),
                # len=0.6,
                # y=0.25,
                thickness=10,
                outlinewidth=0.5,
                outlinecolor="black",
            ),
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Metabolite affinity [mM]",
            zeroline=False,
            type="log",
            dtick="1",
            ticks="inside",
        ),
        yaxis=dict(
            title="Metabolite concentration [mM]",
            zeroline=False,
            type="log",
            dtick="1",
            ticks="inside",
        ),
        height=height * 1.3,
        width=250,
        title="Realizable growth rates",
    )

    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=40,
        buttom_margin=25,
        top_margin=5,
        right_margin=20,
    )
    fig.write_image("plots/simulations/coexistence/fig2e.svg")


def sfig3():
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
        width=200,
        height=200,
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=lm,
        buttom_margin=25,
        top_margin=5,
        right_margin=rm,
        marker_size=4,
    )
    fig.write_image("plots/experiments/sfig3.svg")
