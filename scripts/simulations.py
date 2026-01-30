import numpy as np
from scipy.integrate import odeint
import pandas as pd
from style import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from models import (
    competition as cp,
    niche_creation as nc,
    plot_competition as plot_comp,
    plot_mutual_cf as plot_mutual_cf,
)

# Margins for plots
lm = 10
bm = 10
tm = 10
rm = 10
font_size = 8
line_thickness = 1.2
xs = np.linspace(0, 5000, 5000 * 6)


def parse_params():
    df = dict(pd.read_csv("parameters.csv"))
    params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
    p = params
    return p


def chemostat_acetate_concentration():
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
    fig.write_image("plots/simulations/coexistence/chemostat_acetate_concentration.svg")


def achievable_growth_rate():
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
    fig.write_image("plots/simulations/coexistence/achievable_growth_rate.svg")


def missing_growth_rate():
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
    fig.add_vline(x=0.159)
    fig.write_image("plots/simulations/coexistence/missing_growth_rate.svg")


def metabolite_affinity():
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
    fig.add_trace(
        go.Contour(
            z=zs,
            x=Kms,
            y=Rs,
            showscale=False,
            contours=dict(start=0.07, end=0.07, coloring="none"),
            line=dict(color="black"),
            name="<i>Ct</i>",
            showlegend=False,
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
    fig.write_image("plots/simulations/coexistence/metabolite_affinity.svg")


def km_across_substrates():
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
    fig.write_image("plots/experiments/km_across_substrates.svg")


def simulate_chemostat_community_experiments():
    fig = make_subplots(
        rows=1,
        cols=2,
        horizontal_spacing=0.05,
        subplot_titles=["Cross-feeding", "No cross-feeding"],
        shared_yaxes=True,
    )
    for trace in plot_mutual_cf().data:
        fig.add_trace(trace, row=1, col=1)
    for trace in plot_comp().data:
        fig.add_trace(trace, row=1, col=2)
    fig.update_layout(
        width=190,
        height=180,
        yaxis=dict(title="OD", ticks="inside"),
        xaxis=dict(title="Time [h]", ticks="inside"),
        showlegend=False,
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=20,
        buttom_margin=30,
        top_margin=20,
        right_margin=0,
    )
    fig.write_image(
        "plots/simulations/coexistence/simulate_chemostat_community_experiments.svg"
    )


def simulate_cross_feeding_batch():
    fig = make_subplots(
        rows=1,
        cols=2,
        horizontal_spacing=0.05,
        column_titles=["Ct", "Oa", "Co-culture"],
        shared_yaxes=True,
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
        col=1,
    )

    p = parse_params()
    p["D"] = 0
    # p["N01"] = 0

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
        col=2,
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
        col=2,
    )

    fig.for_each_xaxis(lambda x: x.update(ticks="inside"))
    fig.for_each_yaxis(lambda y: y.update(ticks="inside"))
    fig.update_layout(
        width=190,
        height=180,
        showlegend=False,
        yaxis=dict(title="OD"),
        xaxis2=dict(title="Time [h]"),
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=10,
        buttom_margin=10,
        top_margin=20,
        right_margin=rm,
    )
    fig.write_image("plots/simulations/coexistence/fig3a.svg")

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=M,
            name="M",
            line=dict(color="black", shape="spline"),
        ),
    )
    fig.update_layout(
        width=190,
        height=height,
        showlegend=False,
        xaxis=dict(ticks="inside"),
        yaxis=dict(ticks="inside"),
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=lm,
        buttom_margin=0,
        top_margin=10,
        right_margin=rm,
    )
    fig.write_image(
        "plots/simulations/coexistence/simulate_cross_feeding_batch_metabolite.svg"
    )


def main():
    chemostat_acetate_concentration()
    achievable_growth_rate()
    missing_growth_rate()
    metabolite_affinity()
    km_across_substrates()
    simulate_chemostat_community_experiments()
    simulate_cross_feeding_batch()


main()
