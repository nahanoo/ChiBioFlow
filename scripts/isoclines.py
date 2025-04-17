import numpy as np
import pandas as pd
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots


df = dict(pd.read_csv("parameters.csv"))
params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
p = params

lm = 10
bm = 10
tm = 10
rm = 10
font_size = 8
line_thickness = 1.2
Rs = np.linspace(1e-10, 0.5, 100)


def competition():
    # Competition for acetate
    JCt = [p["v1_1"] * R / (R + p["K1_1"]) for R in Rs]
    JOa = [p["v2_1"] * R / (R + p["K2_1"]) for R in Rs]
    Ds = np.linspace(0, max(JCt), 100)
    R_grid, D_grid = np.meshgrid(Rs, Ds)
    JCt_grid = p["v1_1"] * R_grid / (R_grid + p["K1_1"])
    JOa_grid = p["v2_1"] * R_grid / (R_grid + p["K2_1"])
    R_ratio = np.log(JOa_grid / JCt_grid)
    custom_colorscale = [
        [0, "lightblue"],  # Min value -> Blue
        [0.5, "white"],  # Mid value (zero) -> White
        [1, "lightsalmon"],  # Max value -> Red
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=R_ratio,
            x=Rs,
            y=Ds,
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
                y=0.4,
            ),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=Rs, y=JCt, marker=dict(color=colors["ct"]), name="<i>Ct</i>", opacity=1
        )
    )
    fig.add_trace(
        go.Scatter(x=Rs, y=JOa, marker=dict(color=colors["oa"]), name="<i>Oa</i>")
    )
    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="J [1/h]")
    fig.update_layout(height=height, width=width)
    fig.update_layout(legend_title_text="Isocline of<br>growth rate")
    fig.update_layout(
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
    fig.write_image("plots/isoclines/competition.svg")


def mutual_cf():
    Ts = np.linspace(1e-10, 30, 100)
    R_grid, T_grid = np.meshgrid(Rs, Ts)
    JCt = p["v1_1"] * R_grid / (R_grid + p["K1_1"])
    JOa = p["v2_1"] * R_grid / (R_grid + p["K2_1"]) * T_grid / (T_grid + p["K2_3"])
    J_ratio = np.log10(JOa / JCt)
    custom_colorscale = [
        [0, "#c0bedc"],
        [0.5, "white"],
        [1, "#ecaf80"],
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=J_ratio,
            x=Rs,
            y=Ts,
            colorscale=custom_colorscale,
            zmid=0,
            zmin=-0.1,
            zmax=0.1,
            ncontours=50,
            contours=dict(showlines=False),
            colorbar=dict(
                len=0.8,
                y=0.25,
                thickness=10,
            ),
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
            showlegend=False,
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
            showlegend=False,
        )
    )

    fig.update_layout(
        height=110,
        width=150,
        # legend_title_text="Isocline of<br>growth rate",
        title="Growth rate ratio of Oa to Ct",
        xaxis=dict(
            showgrid=True,
        ),
        yaxis=dict(showgrid=True),
    )
    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="Thiamine [nM]")
    fig = style_plot(
        fig,
        line_thickness=1,
        font_size=font_size,
        left_margin=25,
        buttom_margin=30,
        top_margin=15,
        right_margin=0,
    )
    fig.write_image("plots/isoclines/mutual_cf.svg")


mutual_cf()


def niche_creation():
    Rs = np.linspace(1e-10, 0.5, 1000)
    Ms = np.linspace(1e-10, 0.1, 1000)
    R_grid, M_grid = np.meshgrid(Rs, Ms)
    JCt = p["v1_1"] * R_grid / (R_grid + p["K1_1"]) + p["v1_2"] * M_grid / (
        M_grid + p["K1_2"]
    )
    JOa = p["v2_1"] * R_grid / (R_grid + p["K2_1"])

    J_ratio = np.log10(JOa / JCt)
    custom_colorscale = [
        [0, "lightblue"],  # Min value -> Blue
        [0.5, "white"],  # Mid value (zero) -> White
        [1, "lightsalmon"],  # Max value -> Red
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            z=J_ratio,
            x=Rs,
            y=Ms,
            colorscale=custom_colorscale,
            zmid=0,
            zmin=-0.1,
            zmax=0.1,
            ncontours=50,
            contours=dict(showlines=False),
            colorbar=dict(
                title=dict(
                    text="log<sub>10</sub> ( J<sub><i>Oa</i></sub> / J<sub><i>Ct</i></sub> )",
                    side="right",
                ),
                len=0.5,
                y=0.4,
            ),
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
    fig.update_layout(
        height=height,
        width=width,
        legend_title_text="Isocline of<br>growth rate",
        xaxis=dict(showgrid=False, ticks="outside"),
        yaxis=dict(showgrid=False, ticks="outside"),
    )
    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="Metabolite [mM]")
    fig = style_plot(
        fig,
        line_thickness=line_thickness,
        font_size=font_size,
        left_margin=lm,
        buttom_margin=bm,
        top_margin=tm,
        right_margin=rm,
    )
    fig.write_image("plots/isoclines/niche_creation.svg")


def niche_creation_cf():
    Rs = np.linspace(1e-10, 0.5, 10)
    Ms = np.linspace(1e-10, 0.5, 10)
    Ts = np.linspace(1e-10, p["M3"] / 20, 100)
    R_grid, M_grid, T_grid = np.meshgrid(Rs, Ms, Ts)
    JCt = p["v1_1"] * R_grid / (R_grid + p["K1_1"]) + p["v1_2"] * M_grid / (
        M_grid + p["K1_2"]
    )
    JOa = p["v2_1"] * R_grid / (R_grid + p["K2_1"]) * T_grid / (T_grid + p["K2_3"])
    J_ratio = JOa / JCt
    custom_colorscale = [
        [0, "lightblue"],  # Min value -> Blue
        [0.5, "white"],  # Mid value (zero) -> White
        [1, "lightsalmon"],  # Max value -> Red
    ]
    fig = go.Figure()
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
