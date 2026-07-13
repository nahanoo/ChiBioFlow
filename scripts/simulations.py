import numpy as np
from scipy.integrate import odeint
import pandas as pd
from style import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from chibio_parser import fluorescence_paresr, calibration_csv, cfu_parser
from models import (
    competition as cp,
    niche_creation as nc,
    thiamine_supply as ts,
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
                line=dict(
                    color="black",
                    width=1.5,
                ),
            ),
            dict(
                type="line",
                yref="paper",
                x0=0.15,
                y0=0,
                x1=0.15,
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
        line_thickness=2.5,
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
        yaxis=dict(title="J 1/h", ticks="inside", range=[0, 0.15]),
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
        line_thickness=2.5,
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
            contours=dict(start=0.1, end=0.1, coloring="none"),
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
        height=height * 1.15,
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
        line_thickness=2.5,
    )
    fig.write_image(
        "plots/simulations/coexistence/simulate_chemostat_community_experiments.svg"
    )


def simulate_cross_feeding_batch():
    fig = make_subplots(
        rows=1,
        cols=2,
        horizontal_spacing=0.05,
        column_titles=["Mono-culture", "Co-culture"],
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
    M_prod = Y[:, 3]
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
        line_thickness=2.5,
    )
    fig.write_image("plots/simulations/coexistence/masked_interactions_batch.svg")

    fig = go.Figure()
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=M_prod,
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


def simulate_oa_thiamine_gradient():
    thiamine_colors = {
        0: "#1f77b4",
        0.01: "#ff7f0e",
        0.1: "#2ca02c",
        1: "#d62728",
        10: "#9467bd",
        100: "#8c564b",
        1000: "#e377c2",
        10000: "#7f7f7f",
    }
    p = parse_params()
    p["N02"] = 0.2
    xs = np.linspace(0, 72, 2000)

    fig = go.Figure()
    for M3, color in reversed(thiamine_colors.items()):
        p["M3"] = M3
        Y = odeint(ts, [0, p["N02"], p["M1"], M3], xs, args=(p,))
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=Y[:, 1],
                mode="lines",
                name=f"{M3} nM",
                line=dict(color=color),
                showlegend=True,
            )
        )
    fig.update_layout(
        xaxis=dict(title="Time [h]", ticks="inside", dtick=12),
        yaxis=dict(title="OD", ticks="inside", range=[0, 0.5], dtick=0.1),
        legend=dict(title="Thiamine", font=dict(size=9)),
        width=300,
        height=260,
        title="Oa thiamine gradient simulation",
        showlegend=True,
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=30,
        right_margin=80,
        buttom_margin=30,
        top_margin=35,
    )
    fig.write_image("plots/simulations/dynamics/oa_thiamine_gradient.svg")


def simulate_thiamine_carryover():
    # 8 mL preculture (10000 nM thiamine) + 17 mL thiamine-free medium = 25 mL
    T0 = (8 / 25) * 10000  # = 3200 nM initial thiamine
    N0 = 0.1  # starting OD after dilution
    p = parse_params()
    p["D"] = 0.15
    p["M3"] = 0  # thiamine-free feed

    xs = np.linspace(0, 72, 2000)
    Y = odeint(ts, [0, N0, p["M1"], T0], xs, args=(p,))

    # load experimental data (M3 reactor from 250310_oa_mono)
    Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/250310_oa_mono")
    df_exp = calibration_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250310_oa_mono/calibration.csv", Ms
    )
    df_exp = df_exp[df_exp["reactor"] == "M3"]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df_exp["exp_time"][4:-10],
            y=df_exp["od_calibrated"][4:-10],
            mode="markers",
            name="Experiment",
            marker=dict(color=colors["oa"], size=3),
            showlegend=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 1],
            name="Simulation",
            line=dict(color=colors["oa"]),
            showlegend=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 3],
            name="Thiamine [nM]",
            line=dict(color=colors["blue"], dash="dash"),
            showlegend=True,
            yaxis="y2",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]", ticks="inside", dtick=12),
        yaxis=dict(title="OD", ticks="inside", range=[0, 0.5], dtick=0.1),
        yaxis2=dict(
            title="Thiamine [nM]",
            overlaying="y",
            side="right",
            ticks="inside",
            showgrid=False,
        ),
        legend=dict(font=dict(size=9)),
        width=320,
        height=220,
        title="Oa carryover thiamine chemostat",
        showlegend=True,
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=35,
        right_margin=80,
        buttom_margin=30,
        top_margin=35,
    )
    fig.write_image("plots/simulations/dynamics/oa_thiamine_carryover.svg")
    print(
        f"T0 = {T0:.0f} nM | Final OD: {max(Y[-1, 1], 0):.4f} | Final T: {Y[-1, 3]:.4f} nM"
    )


def thiamine_feed_contamination():
    p = parse_params()
    D = 0.15
    e = "/home/eric/ChiBioFlow/data/260629_oa_washout_ct_oa_no_cs"

    # Load experimental data
    od_data = pd.read_excel(f"{e}/od.ods", engine="odf")
    od_data = od_data[od_data["reactor"] == "M0"].reset_index(drop=True)
    cfus = cfu_parser(e)[0]
    cfus = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "oa")].reset_index(
        drop=True
    )

    # Initial conditions and CFU/OD conversion from t=0
    N0_od = float(od_data.loc[od_data["sample_time"] == 0, "OD"].values[0])
    N0_cfu = float(cfus.loc[cfus["sample_time"] == 0, "average"].values[0])
    cfu_per_od = N0_cfu / N0_od

    # Observed SS OD (last data point) → required feed thiamine T_in
    N_ss_od = float(od_data.iloc[-1]["OD"])
    T_ss = D * p["K2_3"] / (p["v2_1"] - D)
    T_in = T_ss + N_ss_od / p["q2_3"]
    # T_in = 0.5

    print(f"N0 = {N0_od:.3f} OD  |  {N0_cfu:.2e} CFU/mL")
    print(f"cfu_per_od              = {cfu_per_od:.2e}")
    print(f"N_ss = {N_ss_od:.3f} OD")
    print(f"T* (reactor SS)         = {T_ss:.2f} nM")
    print(f"T_in (feed contamination) = {T_in:.2f} nM")

    def ode(y, t, D_val, T_in_val):
        Oa, T = y
        JOa = p["v2_1"] * T / (T + p["K2_3"])
        dOa = JOa * Oa - D_val * Oa
        dT = -JOa * Oa / p["q2_3"] - D_val * T + T_in_val * D_val
        return dOa, dT

    t_carryover = 24  # hours before actual washout begins

    # Phase 0: carryover (flat at N0_od) 0-24 h
    t0 = np.linspace(0, t_carryover, 500)
    N0_flat = np.full(len(t0), N0_od)

    # Phase 1: washout 24-96 h (72 h of dilution with trace T_in)
    t1_rel = np.linspace(0, 96 - t_carryover, 2000)
    Y1 = odeint(ode, [N0_od, T_in], t1_rel, args=(D, T_in))
    t1 = t1_rel + t_carryover

    # Phase 2: dilution stopped at t=96 h
    t2_rel = np.linspace(0, 200, 2000)
    Y2 = odeint(ode, Y1[-1], t2_rel, args=(0.0, 0.0))
    t2 = t2_rel + 96

    fold = Y2[-1, 0] / Y1[-1, 0]
    print(f"CFUs at t=96 h:                       {Y1[-1, 0] * cfu_per_od:.2e}")
    print(f"CFUs after dilution stops (plateau):  {Y2[-1, 0] * cfu_per_od:.2e}")
    print(f"Fold increase after stopping dilution: {fold:.1f}x")

    # --- OD plot ---
    fig_od = go.Figure()
    fig_od.add_trace(
        go.Scatter(
            x=t0,
            y=N0_flat,
            name="Carryover",
            showlegend=True,
            line=dict(color=colors["oa"], dash="dot"),
        )
    )
    fig_od.add_trace(
        go.Scatter(
            x=t1,
            y=Y1[:, 0],
            name="Chemostat",
            showlegend=True,
            line=dict(color=colors["oa"]),
        )
    )
    fig_od.add_trace(
        go.Scatter(
            x=t2,
            y=Y2[:, 0],
            name="Batch",
            showlegend=True,
            line=dict(color=colors["oa"], dash="dash"),
        )
    )
    fig_od.add_trace(
        go.Scatter(
            x=od_data["sample_time"],
            y=od_data["OD"],
            name="Data",
            showlegend=True,
            mode="markers",
            marker=dict(color=colors["oa"], symbol="circle-open", size=8),
        )
    )
    fig_od.add_vline(x=t_carryover, line=dict(color="black", dash="dot", width=1))
    fig_od.add_vline(x=96, line=dict(color="black", dash="dot", width=1))
    fig_od.update_layout(
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD", ticks="inside"),
        width=260,
        height=180,
        title=f"Oa washout (T_in = {T_in:.1f} nM)",
        showlegend=True,
    )
    fig_od = style_plot(
        fig_od,
        font_size=11,
        left_margin=35,
        right_margin=10,
        buttom_margin=30,
        top_margin=30,
    )
    fig_od.write_image("plots/simulations/dynamics/thiamine_feed_contamination_od.svg")

    # --- CFU plot ---
    fig_cfu = go.Figure()
    fig_cfu.add_trace(
        go.Scatter(
            x=t0,
            y=N0_flat * cfu_per_od,
            name="Carryover",
            showlegend=True,
            line=dict(color=colors["oa"], dash="dot"),
        )
    )
    fig_cfu.add_trace(
        go.Scatter(
            x=t1,
            y=Y1[:, 0] * cfu_per_od,
            name="Chemostat",
            showlegend=True,
            line=dict(color=colors["oa"]),
        )
    )
    fig_cfu.add_trace(
        go.Scatter(
            x=t2,
            y=Y2[:, 0] * cfu_per_od,
            name="Batch",
            showlegend=True,
            line=dict(color=colors["oa"], dash="dash"),
        )
    )
    fig_cfu.add_trace(
        go.Scatter(
            x=cfus["sample_time"],
            y=cfus["average"],
            error_y=dict(type="data", array=cfus["stdev"].to_list(), visible=True),
            name="Data",
            showlegend=True,
            mode="markers",
            marker=dict(color=colors["oa"], symbol="circle-open", size=8),
        )
    )
    fig_cfu.add_vline(x=t_carryover, line=dict(color="black", dash="dot", width=1))
    fig_cfu.add_vline(x=96, line=dict(color="black", dash="dot", width=1))
    fig_cfu.update_layout(
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="CFUs/mL", type="log", ticks="inside", exponentformat="power"),
        width=260,
        height=180,
        title=f"Oa washout (T_in = {T_in:.1f} nM)",
        showlegend=True,
    )
    fig_cfu = style_plot(
        fig_cfu,
        font_size=11,
        left_margin=45,
        right_margin=10,
        buttom_margin=30,
        top_margin=30,
    )
    fig_cfu.write_image(
        "plots/simulations/dynamics/thiamine_feed_contamination_cfu.svg"
    )


thiamine_feed_contamination()


def main():
    chemostat_acetate_concentration()
    achievable_growth_rate()
    missing_growth_rate()
    metabolite_affinity()
    km_across_substrates()
    simulate_chemostat_community_experiments()
    simulate_cross_feeding_batch()
