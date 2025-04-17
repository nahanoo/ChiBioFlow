import pandas as pd
import plotly.graph_objects as go
from style import *
import numpy as np
from scipy.stats import linregress


colors = {
    "0.045 OD of Ct": "#440154",  # Dark Purple
    "0.152 OD of Ct": "#3b528b",  # Blue-Purple
    "0.235 OD of Ct": "#21908d",  # Teal
    "0.299 OD of Ct": "#5dc963",  # Light Green
    "0.37 OD of Ct": "#fde725",  # Yellow
}


df = pd.read_csv("data/metadata.csv")
df = df[df["exp_ID"] == "ct_oa_chemostat_project/_oa_in_spent_media_of_ct"]
data = pd.read_csv("data/measurements.csv")
fig = go.Figure()


for i, (lg, OD) in enumerate(zip(df["linegroup"], df["comments"])):
    x = data[lg + "_time"]
    y = data[lg + "_measurement"]
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name=OD,
            legendgroup=OD,
            showlegend=(i <= 4),
            marker=dict(color=colors[OD]),
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
fig = style_plot(fig, line_thickness=1.7)
fig.write_image("plot.svg")


lgs_10000_nM_thiamine = [
    "ct_oa_chemostat_project/_thiamine_gradient_B2",
    "ct_oa_chemostat_project/_thiamine_gradient_C2",
    "ct_oa_chemostat_project/_thiamine_gradient_D2",
]

lgs_ct_037 = [
    "ct_oa_chemostat_project/_oa_in_spent_media_of_ct_B6",
    "ct_oa_chemostat_project/_oa_in_spent_media_of_ct_C6",
    "ct_oa_chemostat_project/_oa_in_spent_media_of_ct_D6",
]
ys = []
for i, lg in enumerate(lgs_ct_037):
    x = data[lg + "_time"]
    y = data[lg + "_measurement"]
    ys.append(y[:30])

y = np.average(ys, axis=0)
slope, intercept, r_value, p_value, std_err = linregress(x[:30].to_numpy(), np.log(y))

fig = go.Figure()


for i, (lg, OD) in enumerate(zip(df["linegroup"], df["comments"])):
    x = data[lg + "_time"]
    y = data[lg + "_measurement"]
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name=OD,
            legendgroup=OD,
            showlegend=(i <= 4),
            marker=dict(color=colors[OD]),
        )
    )
for i, lg in enumerate(lgs_10000_nM_thiamine):
    x = data[lg + "_time"]
    y = data[lg + "_measurement"]
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            name="0.01 mM Thiamine",
            legendgroup=OD,
            showlegend=(i <= 0),
            marker=dict(color="#2d708e"),
        )
    )
"""fig.add_trace(
    go.Scatter(
        x=x[:30],
        y=[0.05 * np.exp(slope * t) for t in x[:30]],
        name=str(slope)[:4] + " 1/h max. growth rate",
        legendgroup=OD,
        opacity=0.5,
        showlegend=True,
        marker=dict(color="blue"),
    )
)"""
fig.update_layout(
    xaxis=dict(
        range=[0, max(x)], showgrid=True, zeroline=True, dtick=6, title="Time [h]"
    ),
    yaxis=dict(
        range=[0, 0.8], showgrid=True, zeroline=True, dtick=0.1, title="Corrected OD"
    ),
)


fig = style_plot(fig, line_thickness=1.7)
fig.write_image("plot_reference.svg")
