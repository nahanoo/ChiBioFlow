import pandas as pd
import json
import plotly.graph_objects as go
import numpy as np
import plotly.express as px
from style import *

colorscale = px.colors.sequential.Viridis  #
n = 4
sampled_colors = [
    colorscale[int(i * (len(colorscale) - 1) / (n - 1))] for i in range(n)
]

f = "master_data.json"
with open(f, "r") as file:
    data = json.load(file)

df = pd.DataFrame(data["P1"])
ct_oa = ["B2", "B3", "B4"]
com_no_cs = ["C2", "C3", "C4"]
com_gluc_c1 = ["D2", "D3", "D4"]
com_gluc_c2 = ["E2", "E3", "E4"]
names = [
    "Ct Oa no CS",
    "Community no CS",
    "Community Glucose C1",
    "Community Glucose C2",
]

fig = go.Figure()
for i, (name, wells) in enumerate(
    zip(names, [ct_oa, com_no_cs, com_gluc_c1, com_gluc_c2])
):
    for well in wells:
        x = np.array(df.at["timepoint", well]) * 10 / 60
        y = np.array(df.at["values", well])
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name=name,
                legendgroup=name,
                showlegend=False if well != wells[0] else True,
                marker=dict(color=sampled_colors[i]),
            )
        )
fig.update_layout(
    xaxis=dict(title="Time [h]", ticks="inside"),
    yaxis=dict(title="OD600", ticks="inside"),
)
fig = style_plot(fig, line_thickness=1.5)
fig.write_image("plot.svg")
fig.show()
