import pandas as pd
import plotly.graph_objects as go
from style import *
from models import ct_mono, oa_mono, parse_params
import numpy as np
from scipy.integrate import odeint
import plotly.io as pio
from scipy.stats import linregress

pio.kaleido.scope.mathjax = None


colors = {
    "0 nM thiamine": "#1f77b4",
    "0.01 nM thiamine": "#ff7f0e",
    "0.1 nM thiamine": "#2ca02c",
    "1 nM thiamine": "#d62728",
    "10 nM thiamine": "#9467bd",
    "100 nM thiamine": "#8c564b",
    "1000 nM thiamine": "#e377c2",
    "10000 nM thiamine": "#7f7f7f",
    "ct": "#7570B3",
    "oa": "#D95F02",
}


def plot_gradient():
    df = pd.read_csv("data/metadata.csv")
    df = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
    ]
    data = pd.read_csv("data/measurements.csv")
    fig = go.Figure()

    for i, (lg, t_conc) in enumerate(zip(df["linegroup"], df["comments"])):
        x = data[lg + "_time"]
        y = data[lg + "_measurement"]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name=" ".join(t_conc.split(" ")[:2]),
                legendgroup=t_conc,
                showlegend=(i <= 7),
                marker=dict(color=colors[t_conc]),
            )
        )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)], showgrid=True, zeroline=True, dtick=6, title="Time [h]"
        ),
        yaxis=dict(
            range=[0, 0.8],
            showgrid=True,
            zeroline=True,
            dtick=0.1,
            title="Corrected OD",
        ),
        title="Thiamine gradient <i>O. anthropi</i>",
    )
    fig = style_plot(fig, line_thickness=1.7)
    fig.write_image("thiamine_gradient.pdf")


df = pd.read_csv("data/metadata.csv")
df_ct = df[
    (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
    & (df["species"] == "Comamonas testosteroni")
    & (df["comments"] == "10000 nM thiamine")
]
df_oa = df[
    (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
    & (df["species"] == "Ochrobactrum anthropi")
    & (df["comments"] == "10000 nM thiamine")
]
data = pd.read_csv("data/measurements.csv")
p = parse_params()
p["D"] = 0
xs = np.linspace(0, 72, 200)
Y_ct = odeint(ct_mono, [p["N01"], p["M1"]], xs, args=(p,))
Y_oa = odeint(oa_mono, [p["N02"], p["M1"]], xs, args=(p,))

fig = go.Figure()
for i, lg in enumerate(df_ct["linegroup"]):
    x = data[lg + "_time"]
    y = data[lg + "_measurement"]
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="<i>C. testosteroni</i>",
            legendgroup="<i>C. testosteroni</i>",
            showlegend=(i == 0),
            marker=dict(color=colors["ct"]),
        )
    )
fig.add_trace(
    go.Scatter(
        x=xs,
        y=Y_ct[:, 0],
        name="<i>C. testosteroni</i> model",
        marker=dict(color=colors["ct"]),
        line=dict(dash="dash"),
    )
)
for i, lg in enumerate(df_oa["linegroup"]):
    x = data[lg + "_time"]
    y = data[lg + "_measurement"]
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="<i>O. anthropi</i>",
            legendgroup="<i>O. anthropi</i>",
            showlegend=(i == 0),
            marker=dict(color=colors["oa"]),
        )
    )

fig.add_trace(
    go.Scatter(
        x=xs,
        y=Y_oa[:, 0],
        name="<i>O. anthropi</i> model",
        marker=dict(color=colors["oa"]),
        line=dict(dash="dash"),
    )
)

y_oa = [data[lg + "_measurement"] for lg in df_ct["linegroup"]]
y_oa = np.array(y_oa).mean(axis=0)[:36]
x_oa = data[df_oa["linegroup"].iloc[0] + "_time"][:36]
slope, intercept, r_value, p_value, std_err = linregress(x_oa, np.log(y_oa))
fit = [np.exp(slope * i + intercept) for i in x_oa]
fig.add_trace(
    go.Scatter(
        x=x_oa,
        y=fit,
        name="Fit",
        marker=dict(color="black"),
        line=dict(dash="dash"),
    )
)


fig.update_layout(
    xaxis=dict(
        range=[0, max(x)], showgrid=True, zeroline=True, dtick=6, title="Time [h]"
    ),
    yaxis=dict(
        range=[0, 0.8], showgrid=True, zeroline=True, dtick=0.1, title="Corrected OD"
    ),
)

fig.update_layout(width=width, height=height)
fig = style_plot(fig, line_thickness=1.7)
fig.write_image("ct_oa_mono_culture.svg")
