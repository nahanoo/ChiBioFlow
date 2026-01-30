import plotly.graph_objects as go
from chibio_parser import cfu_parser
from plotly.subplots import make_subplots
from style import *
import pandas as pd

cfus = cfu_parser("/home/eric/ChiBioFlow/data/250822_chemostats")[0]


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
    "ms": "Ml",
    "at": "At",
}

od = pd.read_excel("od_readings/od.xlsx")


def plot_experiment(cfus, od):
    fig = make_subplots(
        rows=4,
        cols=2,
        shared_xaxes=True,
        shared_yaxes=True,
        subplot_titles=["No CS", "No CS", "", "", "Glucose C1", "Glucose C2"],
        vertical_spacing=0.05,
        row_heights=[0.3, 0.2, 0.3, 0.2],
    )
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
                    error_y=dict(
                        type="data", array=df["stdev"].to_list(), visible=True
                    ),
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
        range=[0, 0.2],
        title="OD",
        row=2,
        col=1,
    )
    fig.update_yaxes(
        range=[0, 0.2],
        title="OD",
        row=2,
        col=2,
    )
    fig.update_yaxes(
        range=[0, 0.7],
        title="OD",
        row=4,
        col=1,
    )
    fig.update_yaxes(
        range=[0, 0.7],
        title="OD",
        row=4,
        col=2,
    )
    fig.update_xaxes(title=None)
    fig.update_xaxes(title="Time [h]", row=4, col=1)
    fig.update_xaxes(title="Time [h]", row=4, col=2)

    fig.update_layout(height=800, width=1000)
    fig = style_plot(fig, marker_size=10, font_size=11, line_thickness=2)
    fig.write_image("plot.svg")


plot_experiment(cfus, od)


def snorre(cfus, od):
    fig = go.Figure()
    cfus = cfus[cfus["reactor"] == "M1"]
    for s in cfus["species"].unique():
        df = cfus[cfus["species"] == s]
        fig.add_trace(
            go.Scatter(
                x=df["sample_time"],
                y=df["average"],
                mode="lines+markers",
                name=names[s],
                legendgroup=s,
                marker=dict(color=colors[s]),
                error_y=dict(type="data", array=df["stdev"].to_list(), visible=True),
            )
        )
    fig.update_yaxes(
        type="log",
        range=[2, 11],
        dtick=1,
        title="CFUs/mL",
        tickformat=".0e",
    )
    fig.update_xaxes(title="Time [h]")
    fig.update_layout(height=350, width=500, title="No carbon source")
    fig = style_plot(fig, marker_size=10, font_size=11, line_thickness=2)
    fig.write_image("community_no_cs.svg")
    cfus = cfus[
        ["species", "reactor", "sample_time", "dilution", "count", "average", "stdev"]
    ]
    cfus.to_csv("community_no_cs.csv", index=False)
    od.index = od["time"]
    od = od[["M1"]]
    od.to_csv("od_commmunity_no_cs.csv")


legend = {"M1": "No CS", "M2": "Glucose", "M3": "Second reactor<br>in chain"}
colors = ["#7570b3", "#1b9e77"]
df = pd.DataFrame(columns=["time", "species", "cfus", "reactor"])
for i, r in enumerate(["M2", "M3"]):
    for j, s in enumerate(["ct", "oa", "ms", "at"]):
        for t in cfus["sample_time"].unique():
            if t == 0.0:
                continue
            avg = cfus[
                (cfus["reactor"] == r)
                & (cfus["species"] == s)
                & (cfus["sample_time"] == t)
            ]["average"].values[0]
            df.loc[len(df)] = [t, s, avg, r]

fig = go.Figure()
for i, r in enumerate(df["reactor"].unique()):
    for s in df["species"].unique():
        dff = df[(df["species"] == s) & (df["reactor"] == r)]
        y = dff[(dff["species"] == s) & (dff["reactor"] == r)]["cfus"].to_numpy()
        x = [names[s]] * len(y)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                marker=dict(color=colors[i]),
                mode="markers",
                offsetgroup=r,
                name=legend[r],
                legendgroup=r,
                showlegend=(True if s == "ct" else False),
            )
        )
fig.update_yaxes(
    type="log",
    range=[5, 10],
    dtick=1,
    title="CFUs/mL",
    tickformat=".0e",
)
fig.update_xaxes(title="Species")
fig.update_layout(
    height=400,
    width=400,
    scattermode="group",
    boxgap=0.2,
)
fig = style_plot(
    fig,
    marker_size=10,
    font_size=11,
)
fig.write_image("chain_comparison.svg")
