import numpy as np
from scipy.integrate import odeint
import pandas as pd
import plotly.graph_objects as go
from style import *

df = dict(pd.read_csv("parameters.csv"))
params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
p = params

width, height = 300, 250


def competition():
    # Competition for acetate
    Rs = np.linspace(0, p["M1"], 1000)
    JCt = [p["v1_1"] * R / (R + p["K1_1"]) for R in Rs]
    JOa = [p["v2_1"] * R / (R + p["K2_1"]) for R in Rs]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=Rs, y=JCt, marker=dict(color=colors["ct"]), name="<i>Ct</i>")
    )
    fig.add_trace(
        go.Scatter(x=Rs, y=JOa, marker=dict(color=colors["oa"]), name="<i>Oa</i>")
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
    fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="Thiamine [mM]")
    fig = style_plot(fig)
    fig.write_image("plots/isoclines/mutual_cf.svg")


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


from plotly.subplots import make_subplots


def niche_creation_cf():
    pass


Rs = np.linspace(0, 0.5, 100)
Ms = np.linspace(0, 0.5, 100)
Ts = np.linspace(0, p["M3"] / 10, 100)
Ct_R_grid, Ct_M_grid = np.meshgrid(Rs, Ms)
Oa_R_grid, Oa_T_grid = np.meshgrid(Rs, Ts)
JCt = p["v1_1"] * Ct_R_grid / (Ct_R_grid + p["K1_1"]) + p["v1_2"] * Ct_M_grid / (
    Ct_M_grid + p["K1_2"]
)
JOa = (
    p["v2_1"]
    * Oa_R_grid
    / (Oa_R_grid + p["K2_1"])
    * Oa_T_grid
    / (Oa_T_grid + p["K2_3"])
)

fig = make_subplots(specs=[[{"secondary_y": True}]])
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
        y=Ms,
        showscale=False,
        contours=dict(start=0.15, end=0.15, size=0.1, coloring="none"),
        line=dict(color=colors["ct"]),
        name="<i>Ct</i>",
    ),
    secondary_y=True,
)
fig.update_layout(height=height, width=width)
fig.update_xaxes(title="Acetate [mM]"), fig.update_yaxes(title="Metabolite [mM]")
fig = style_plot(fig)
fig.write_image("plots/isoclines/niche_creation_cf.svg")


niche_creation_cf()
