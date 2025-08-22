import plotly.graph_objects as go
from chibio_parser import cfu_parser
from plotly.subplots import make_subplots
from style import *
import pandas as pd

cfus = cfu_parser("/home/eric/ChiBioFlow/data/250623_snorre_batch1")[0]

fig = make_subplots(
    rows=4,
    cols=2,
    shared_xaxes=True,
    shared_yaxes=True,
    subplot_titles=["Ribose", "Acetate", "", "", "Histidine", "Glutaric Acid"],
    vertical_spacing=0.05,
    row_heights=[0.3, 0.2, 0.3, 0.2],
)

colors = {
    "blue": "#000080",
    "ct": "#e1812c",
    "oa": "#c03d3e",
    "ms": "#3274a1",
    "at": "#3b923b",
    "Spent media Ct": "#1B9E77",
    "Spent media Oa": "#E7298A",
    "H20": "gray",
}

names = {
    "ct": "Ct",
    "oa": "Oa",
    "ms": "Ms",
    "at": "At",
}

od = pd.read_excel("od_readings/od.xlsx")

for i, r in enumerate(cfus["reactor"].unique()):
    for j, s in enumerate(cfus["species"].unique()):
        row = 2 * (i // 2) + 1
        col = i % 2 + 1
        df = cfus[(cfus["reactor"] == r) & (cfus["species"] == s)]
        fig.add_trace(
            go.Scatter(
                x=df["sample_time"],
                y=df["average"],
                mode="lines+markers",
                name=names[s],
                legendgroup=s,
                marker=dict(color=colors[s]),
                error_y=dict(type="data", array=df["stdev"].to_list(), visible=True),
                showlegend=(True if i == 0 else False),
            ),
            row=row,
            col=col,
        )

for i, r in enumerate(od.columns[1:]):
    fig.add_trace(
        go.Scatter(
            x=od["time"],
            y=od[r],
            mode="lines+markers",
            name="OD",
            legendgroup="od",
            marker=dict(color=colors["H20"]),
            showlegend=(True if i == 0 else False),
        ),
        row=2 * (i // 2) + 2,
        col=i % 2 + 1,
    )

fig.update_yaxes(
    type="log",
    range=[2, 11],
    dtick=1,
    title="CFUs/mL",
    tickformat=".0e",
    row=1,
    col=1,
)
fig.update_yaxes(
    type="log",
    range=[2, 11],
    dtick=1,
    title="CFUs/mL",
    tickformat=".0e",
    row=1,
    col=2,
)
fig.update_yaxes(
    type="log",
    range=[2, 11],
    dtick=1,
    title="CFUs/mL",
    tickformat=".0e",
    row=3,
    col=1,
)
fig.update_yaxes(
    type="log",
    range=[2, 11],
    dtick=1,
    title="CFUs/mL",
    tickformat=".0e",
    row=3,
    col=2,
)

# OD axes

fig.update_yaxes(
    range=[0, 1.5],
    title="OD",
    row=2,
    col=1,
)
fig.update_yaxes(
    range=[0, 1.5],
    title="OD",
    row=2,
    col=2,
)
fig.update_yaxes(
    range=[0, 0.5],
    title="OD",
    row=4,
    col=1,
)
fig.update_yaxes(
    range=[0, 0.5],
    title="OD",
    row=4,
    col=2,
)
fig.update_xaxes(title=None)
fig.update_xaxes(title='Time [h]',row=4,col=1)
fig.update_xaxes(title='Time [h]',row=4,col=2)

fig.update_layout(height=800, width=1000)
fig = style_plot(fig, marker_size=10, font_size=11, line_thickness=2)
fig.write_image("plot.svg")

reactor_cs = {'M0':'ribose','M1':'acetate','M2':'histidine','M3':'glutaric acid'}
cfus["reactor"] = cfus["reactor"].map(reactor_cs)
cfus = cfus[['species','reactor','sample_time','dilution','count','average','stdev']]
