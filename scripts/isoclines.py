import numpy as np
import pandas as pd
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots


df = dict(pd.read_csv("parameters.csv"))
params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
p = params

width, height = 300, 250


def competition():
    # Competition for acetate
    Rs = np.linspace(0, p["M1"] / 10, 1000)
    JCt = [p["v1_1"] * R / (R + p["K1_1"]) for R in Rs]
    JOa = [p["v2_1"] * R / (R + p["K2_1"]) for R in Rs]
    R_star_1_1 = -p["D"] * p["K1_1"] / (p["D"] - p["v1_1"])
    R_star_2_1 = -p["D"] * p["K2_1"] / (p["D"] - p["v2_1"])
    Rs_2_1 = np.linspace(0, R_star_2_1, 100)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=Rs, y=JCt, marker=dict(color=colors["ct"]), name="<i>Ct</i>")
    )
    fig.add_trace(
        go.Scatter(x=Rs, y=JOa, marker=dict(color=colors["oa"]), name="<i>Oa</i>")
    )
    fig.add_trace(
        go.Scatter(
            x=[0, R_star_1_1],
            y=[p["D"], p["D"]],
            marker=dict(color="black"),
            line=dict(dash="dot"),
            name="R*",
            showlegend=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[R_star_1_1, R_star_1_1],
            y=[0, p["D"]],
            marker=dict(color="black"),
            line=dict(dash="dot"),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[0, R_star_2_1],
            y=[p["D"], p["D"]],
            marker=dict(color="black"),
            line=dict(dash="dot"),
            name="R*",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[R_star_2_1, R_star_2_1],
            y=[0, p["D"]],
            marker=dict(color="black"),
            line=dict(dash="dot"),
            showlegend=False,
        )
    )

    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="J [1/h]")
    fig.update_layout(height=height, width=width)
    fig = style_plot(fig)
    fig.write_image("plots/isoclines/competition.svg")


def mutual_cf():
    Rs = np.linspace(0, 0.5, 1000)
    Ts = np.linspace(0, p["M3"] / 10, 1000)
    R_grid, T_grid = np.meshgrid(Rs, Ts)
    JCt = p["v1_1"] * R_grid / (R_grid + p["K1_1"])
    JOa = p["v2_1"] * R_grid / (R_grid + p["K2_1"]) * T_grid / (T_grid + p["K2_3"])
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=JOa,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, size=0.1, coloring="none"),
            line=dict(color=colors["oa"]),
            name="<i>Oa</i>",
        )
    )
    fig.add_trace(
        go.Contour(
            z=JCt,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, size=0.1, coloring="none"),
            line=dict(color=colors["ct"]),
            name="<i>Ct</i>",
        )
    )
    fig.update_layout(height=height, width=width)
    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="Thiamine [µM]")
    fig = style_plot(fig)
    fig.write_image("plots/isoclines/mutual_cf.svg")
    J_ration = np.log10(JOa / JCt)
    custom_colorscale = [
        [0, "blue"],  # Min value -> Blue
        [0.5, "white"],  # Mid value (zero) -> White
        [1, "red"],  # Max value -> Red
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=J_ration,
            x=Rs,
            y=Ts,
            colorscale=custom_colorscale,
            zmid=0,
            zmin=-0.5,
            zmax=0.5,
            ncontours=50,
            contours=dict(showlines=False),
            colorbar=dict(
                title=dict(
                    text="log<sub>10</sub> ( J<sub><i>Oa</i></sub> / J<sub><i>Ct</i></sub> )",
                    side="right",
                ),
                len=0.5,
            ),
        )
    )
    fig.add_trace(
        go.Contour(
            z=JOa,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, size=0.1, coloring="none"),
            line=dict(color=colors["oa"]),
            name="<i>Oa</i>",
        )
    )
    fig.add_trace(
        go.Contour(
            z=JCt,
            x=Rs,
            y=Ts,
            showscale=False,
            contours=dict(start=0.15, end=0.15, size=0.1, coloring="none"),
            line=dict(color=colors["ct"]),
            name="<i>Ct</i>",
        )
    )
    fig.update_layout(height=height, width=width)
    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="Thiamine [µM]")
    fig = style_plot(fig, font_size=8)
    fig.write_image("plots/isoclines/mutual_cf_ratio.svg")


def niche_creation():
    Rs = np.linspace(0, 0.5, 1000)
    Ms = np.linspace(0, 0.5, 1000)
    R_grid, M_grid = np.meshgrid(Rs, Ms)
    JCt = p["v1_1"] * R_grid / (R_grid + p["K1_1"]) + p["v1_2"] * M_grid / (
        M_grid + p["K1_2"]
    )
    JOa = p["v2_1"] * R_grid / (R_grid + p["K2_1"])

    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=JOa,
            x=Rs,
            y=Ms,
            showscale=False,
            contours=dict(start=0.15, end=0.15, size=0.1, coloring="none"),
            line=dict(color=colors["oa"]),
            name="<i>Oa</i>",
        )
    )
    fig.add_trace(
        go.Contour(
            z=JCt,
            x=Rs,
            y=Ms,
            showscale=False,
            contours=dict(start=0.15, end=0.15, size=0.1, coloring="none"),
            line=dict(color=colors["ct"]),
            name="<i>Ct</i>",
        )
    )
    fig.update_layout(height=height, width=width)
    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="Metabolite [mM]")
    fig = style_plot(fig)
    fig.write_image("plots/isoclines/niche_creation.svg")


def niche_creation_cf():
    Rs = np.linspace(0, 0.5, 10)
    Ms = np.linspace(0, 0.5, 10)
    Ts = np.linspace(0, p["M3"] / 20, 100)
    R_grid, M_grid, T_grid = np.meshgrid(Rs, Ms, Ts)
    JCt = p["v1_1"] * R_grid / (R_grid + p["K1_1"]) + p["v1_2"] * M_grid / (
        M_grid + p["K1_2"]
    )
    JOa = p["v2_1"] * R_grid / (R_grid + p["K2_1"]) * T_grid / (T_grid + p["K2_3"])

    fig = go.Figure()
    fig.add_trace(
        go.Isosurface(
            x=R_grid.flatten(),
            y=M_grid.flatten(),
            z=T_grid.flatten(),
            value=JCt.flatten(),
            isomin=0.15,
            isomax=0.15,
            surface_count=1,
            colorscale=[[0, colors["ct"]], [1, colors["ct"]]],
            name="<i>Ct</i>",
            showscale=False,
            showlegend=True,
            opacity=0.5,
        )
    )
    fig.add_trace(
        go.Isosurface(
            x=R_grid.flatten(),
            y=M_grid.flatten(),
            z=T_grid.flatten(),
            value=JOa.flatten(),
            isomin=0.15,
            isomax=0.15,
            surface_count=1,
            colorscale=[[0, colors["oa"]], [1, colors["oa"]]],
            name="<i>Oa</i>",
            showlegend=True,
            showscale=False,
            opacity=0.5,
        ),
    )
    fig.add_trace(
        go.Scatter3d(
            x=[0.08 + 0.02],
            y=[0.019],
            z=[0.22],
            marker=dict(color="black", opacity=1),
            showlegend=False,
        )
    )
    fig.update_layout(
        scene=dict(
            zaxis=dict(title="Thamine [µM]"),
            xaxis=dict(title="Acetate [mM]"),
            yaxis=dict(title="Metabolite [mM]"),
        )
    )
    fig.update_layout(
        height=height,
        width=width,
        scene_camera=dict(eye=dict(x=2.0, y=2.0, z=2), center=dict(x=0.3, y=0, z=-0.2)),
    )
    fig = style_plot(fig, top_margin=0, left_margin=0, buttom_margin=0, font_size=7)
    fig.write_image("plots/isoclines/niche_creation_cf.svg")
