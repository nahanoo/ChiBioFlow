import pandas as pd
import plotly.graph_objects as go
from style import *
import numpy as np


df = pd.read_excel("250313_oa_chemostat_samples.xlsx")
fig = go.Figure()
fig.add_trace(
    go.Scatter(x=df["Time"], y=df["C2"], name="C1", marker=dict(color="blue"))
)
fig.add_trace(go.Scatter(x=df["Time"], y=df["C3"], name="C2", marker=dict(color="red")))
fig.add_trace(
    go.Scatter(x=df["Time"], y=df["C4"], name="C3", marker=dict(color="green"))
)
fig.add_trace(
    go.Scatter(x=df["Time"], y=df["C5"], name="C4", marker=dict(color="orange"))
)
fig.update_layout(
    title="Chemostat samples",
    xaxis=dict(
        range=[0, max(df["Time"])],
        showgrid=True,
        zeroline=True,
        dtick=6,
        title="Time [h]",
    ),
    yaxis=dict(
        range=[0, 0.3],
        showgrid=True,
        zeroline=True,
        dtick=0.05,
        title="OD",
    ),
    width=width,
    height=height,
)
fig = style_plot(fig, buttom_margin=80, left_margin=60)
fig.write_image("chemostat_samples.svg")
