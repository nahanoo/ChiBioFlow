import numpy as np
import pandas as pd
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots
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
