import numpy as np
import pandas as pd
import plotly.graph_objects as go
from style import *
from joblib import Parallel, delayed


df = dict(pd.read_csv("parameters.csv"))
params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
p = params


def ct_oa_isocline():
    Ts = np.geomspace(1e-2, 1e4)
    Rs = np.geomspace(1e-2, 10)

    R_grid, T_grid = np.meshgrid(Rs, Ts)
    JCt = p["v1_1"] * R_grid / (R_grid + p["K1_1"])
    JOa = p["v2_1"] * R_grid / (R_grid + p["K2_1"]) * T_grid / (T_grid + p["K2_3"])
    JOat = p["v2_1"] * R_grid / (R_grid + p["K2_1"])
    J_ratio = np.log10(JOa / JCt)

    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=J_ratio,
            x=Rs,
            y=Ts,
            colorscale=colors_heatmap,
            zmid=0,
            zmin=-1,
            zmax=1,
            ncontours=50,
            contours=dict(showlines=False),
            colorbar=dict(
                len=0.6,
                y=0.2,
                thickness=10,
                outlinewidth=0.5,
                outlinecolor="black",
            ),
        )
    )
    fig.update_layout(
        autosize=False,
        height=150,
        width=200,
        yaxis=dict(
            title="Thiamine [nM]",
            zeroline=True,
            showgrid=True,
            ticks="inside",
            type="log",
            dtick=1,
            exponentformat="power",
            showexponent="all",
        ),
        xaxis=dict(
            title="Acetate [mM]",
            zeroline=False,
            showgrid=False,
            ticks="inside",
            type="log",
            dtick=1,
        ),
    )
    fig = style_plot(
        fig,
        line_thickness=1.5,
        font_size=11,
        left_margin=30,
        buttom_margin=25,
        top_margin=15,
        right_margin=0,
    )
    fig.add_trace(
        go.Contour(
            z=JOat,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, coloring="none"),
            line=dict(
                color="black",
                width=3.5,
                dash="4px 1px",
            ),
            showlegend=False,
        )
    )

    fig.add_trace(
        go.Contour(
            z=JOat,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, coloring="none"),
            line=dict(color=colors["oa"], width=1.5),
            name="<i>Oa</i>",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Contour(
            z=JCt,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, coloring="none"),
            line=dict(color="black", width=3.5),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Contour(
            z=JCt,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, coloring="none"),
            line=dict(color=colors["ct"], width=1.5),
            name="<i>Ct</i>",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Contour(
            z=JOa,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, coloring="none"),
            line=dict(
                color="black",
                width=3.5,
            ),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Contour(
            z=JOa,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, coloring="none"),
            line=dict(color=colors["oa"], width=1.5),
            name="<i>Oa</i>",
            showlegend=False,
        )
    )
    fig = square_panel_by_height(fig, height_px=150)
    fig.write_image("plots/isoclines/ct_oa_isoclines.svg")


def compute_ratio(D, alpha, p_base, xs, thiamine_supplied=True):
    p = p_base.copy()
    p["D"] = D
    if thiamine_supplied:
        p["M3"] = alpha
        Y = odeint(ts, [p["N01"], p["N02"], p["M1"], p["M3"]], xs, args=(p,))
    else:
        p["q1_3"] = alpha
        Y = odeint(mc, [p["N01"], p["N02"], p["M1"], 0], xs, args=(p,))
    Ct, Oa, R, T = Y[-1]
    if Ct <= 1e-6:
        Ct = 0
    if Oa <= 1e-6:
        Oa = 0
    if (Ct == 0) and (Oa == 0):
        return np.nan, T  # use NaN instead of None for numeric matrix
    else:
        return Oa / (Ct + Oa), T


def coexistence_sweep_thiamine_added():
    p_base = parse_params()
    Ds = np.linspace(0, 0.3, 100)
    alphas = np.linspace(1, 100, 100)

    # Create full parameter grid
    param_grid = [
        (i, j, D, alpha) for i, D in enumerate(Ds) for j, alpha in enumerate(alphas)
    ]

    # Run in parallel
    results = Parallel(n_jobs=-1, verbose=1)(
        delayed(compute_ratio)(D, alpha, p_base, xs) for (_, _, D, alpha) in param_grid
    )

    # Reconstruct result matrix
    zs = np.array([res[0] for res in results]).reshape(len(Ds), len(alphas))

    # Plot
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
                thickness=10,
            ),
            showscale=False,
        )
    )

    fig.update_xaxes(
        title="Thiamine supply concentration [nM]", zeroline=False, ticks="inside"
    )
    fig.update_yaxes(
        title="Dilution rate [1/h]", zeroline=False, showgrid=False, ticks="inside"
    )
    fig.update_layout(height=150, width=170, title="Thiamine supplied")
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=11,
        left_margin=20,
        buttom_margin=25,
        top_margin=20,
        right_margin=10,
    )
    fig.write_image("plots/simulations/coexistence/coexistence_thiamine_supplied.svg")


def coexistence_sweep_thiamine_free():
    p_base = parse_params()
    Ds = np.linspace(0, 0.3, 200)
    alphas = np.linspace(0.0002, 1, 200)

    param_grid = [
        (i, j, D, alpha) for i, D in enumerate(Ds) for j, alpha in enumerate(alphas)
    ]

    results = Parallel(n_jobs=-1, verbose=1)(
        delayed(compute_ratio)(D, alpha, p_base, xs, thiamine_supplied=False)
        for (_, _, D, alpha) in param_grid
    )

    # Reconstruct matrices
    ratios = np.array([r for r, T in results]).reshape(len(Ds), len(alphas))
    Ts = np.array([T for r, T in results]).reshape(len(Ds), len(alphas))

    fig = go.Figure()

    fig.add_trace(
        go.Contour(
            z=ratios,
            x=Ts[0],  # assumes all rows of Ts have the same alpha ordering
            y=Ds,
            colorscale=colors_heatmap,
            ncontours=50,
            zmid=0.5,
            zmin=0,
            zmax=1,
            contours=dict(showlines=False),
            colorbar=dict(
                title=dict(text="<i>Oa</i> fraction", side="right", font=dict(size=8)),
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

    fig.write_image("plots/simulations/coexistence/coexistence_cross_feeding.svg")
