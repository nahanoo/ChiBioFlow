from chibio_parser import cfu_parser
from style import style_plot
import plotly.graph_objects as go
import pandas as pd

colors = {
    "blue": "#000080",
    "ct": "#e1812c",
    "oa": "#c03d3e",
    "ml": "#3274a1",
    "at": "#3b923b",
    "Spent media Ct": "#1B9E77",
    "Spent media Oa": "#E7298A",
    "H20": "gray",
}


def plot_reactor(cfus, fname, title):
    cfus.replace({"species": {"ms": "ml"}}, inplace=True)
    fig = go.Figure()
    for s, df in cfus.groupby("species"):
        fig.add_trace(
            go.Scatter(
                x=df["sample_time"],
                y=df["average"],
                mode="lines+markers",
                name=s[0].upper() + s[1:],
                marker=dict(color=colors[s]),
                error_y=dict(type="data", array=df["stdev"].to_list(), visible=True),
                showlegend=True,
            )
        )
    fig.update_layout(
        xaxis=dict(title="Time (hours)", ticks="inside"),
        yaxis=dict(title="CFU/mL", type="log", range=[2, 11], ticks="inside"),
        title=title,
        legend=dict(title="Species", x=0.01, y=0.99),
    )
    fig = style_plot(fig, line_thickness=2.5, marker_size=8, font_size=12)
    fig.write_image(fname)


def no_carbon_source():
    e = "250822_chemostats"
    cfus, o = cfu_parser(e)
    cfus = cfus[cfus["reactor"] == "M1"]
    cfus["reactor"] = "no carbon source"
    cfus = cfus[
        ["species", "reactor", "sample_time", "average", "stdev", "count", "dilution"]
    ]
    plot_reactor(cfus, "plots/cfus_no_cs.svg", "No carbon source")
    cfus.to_csv("dataframes/cfus_no_cs.csv", index=False)
    df = pd.read_excel("../250822_chemostats/od_readings/od.xlsx")
    df = df[["time", "M1"]]
    df.columns = ["sample_time", "od"]
    df.to_csv("dataframes/od_no_cs.csv", index=False)


cfus = cfu_parser("250623_snorre_batch1")[0]
reactor_to_cs = {
    "M0": "ribose",
    "M1": "acetate",
    "M2": "histidine",
    "M3": "glutaric acid",
}
cfus.replace({"reactor": reactor_to_cs}, inplace=True)
cfus.replace({"species": {"ms": "ml"}}, inplace=True)
for medium, df in cfus.groupby("reactor"):
    fname = f"plots/cfus_{medium.replace(' ', '_')}.svg"
    title = medium
    # plot_reactor(df, fname, title)
    df = df[
        ["species", "reactor", "sample_time", "average", "stdev", "count", "dilution"]
    ]
    csv_fname = f"dataframes/cfus_{medium.replace(' ', '_')}.csv"
    # df.to_csv(csv_fname, index=False)
df = pd.read_excel("../250623_snorre_batch1/od_readings/od.xlsx")
df.columns = ["sample_time", "ribose", "acetate", "histidine", "glutaric acid"]
for medium in df.columns[1:]:
    csv_fname = f"dataframes/od_{medium.replace(' ', '_')}.csv"
    df_medium = df[["sample_time", medium]]
    df_medium.columns = ["sample_time", "od"]
    df_medium.to_csv(csv_fname, index=False)
