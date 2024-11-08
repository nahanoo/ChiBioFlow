import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots
import pandas as pd
from os.path import join, split
from os import listdir
import numpy as np
from glob import glob

colors = {
    "ct": "#7570B3",
    "oa": "#D95F02",
    "<i>A. tumefaciens</i>": "#1B9E77",
    "<i>O. anthropi</i>": "#D95F02",
}

abb = {
    "at": "<i>A. tumefaciens</i>",
    "oa": "<i>O. anthropi</i>",
    "<i>A. tumefaciens</i>": "at",
    "<i>O. anthropi</i>": "oa",
}


dirs = [
    "240818_ct_oa_thiamine",
    "240818_ct_oa_thiamine_restart",
    "240818_ct_oa_thiamine_restart_2",
]


def style_plot(fig, marker_size=3, top_margin=20, font_size=14):
    """Style function for figures setting fot size and true black color."""
    fig.update_layout(
        {
            "plot_bgcolor": "rgb(168, 168, 168)",
            "paper_bgcolor": "rgb(207, 207, 207)",
        }
    )
    for d in fig["data"]:
        d["marker"]["size"] = marker_size
        d["line"]["width"] = marker_size
    # Font size
    j = font_size
    fig.update_layout(font={"size": j, "color": "black"})
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = j
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = j
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = j
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"
    fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=j, color="black")))
    fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=j, color="black")))
    fig.for_each_yaxis(lambda axis: axis.update(gridcolor="black"))
    fig.for_each_xaxis(lambda axis: axis.update(gridcolor="black"))

    fig.update_layout(
        margin=dict(l=60, r=20, t=top_margin, b=20),
        hoverlabel=dict(font_size=font_size),
    )
    fig.update_yaxes(title_standoff=10)
    fig.update_xaxes(title_standoff=10)
    return fig


def parse_data():
    dir = join("/", "home", "eric", "ChiBioFlow", "data")
    df = fluorescence_paresr(
        join(dir, "at_oa/240921_ct_mono_batch"),
        down_sample=False,
        window_size=20,
        od_filter=1,
    )

    M0 = df[df["reactor"] == "M0"]
    M1 = df[df["reactor"] == "M1"]
    M2 = df[df["reactor"] == "M2"]
    abs = pd.DataFrame(columns=["exp_time", "average", "stdev", "read", "species"])
    fs = listdir(join(dir, "at_oa/240921_ct_mono_batch", "plate_reader", "Absorbance"))
    for f in sorted(fs):
        t = float(f[:-5])
        data = pd.read_excel(
            join(
                dir,
                "at_oa/240921_ct_mono_batch",
                "plate_reader",
                "Absorbance",
                f,
            ),
            index_col=0,
        )
        avg = np.average(data.loc["A", [1, 2, 3]])
        stdev = np.std(data.loc["A", [1, 2, 3]])
        abs.loc[len(abs)] = [t, avg, stdev, "Absorbance", "Community 1"]

        avg = np.average(data.loc["A", [4, 5, 6]])
        stdev = np.std(data.loc["A", [4, 5, 6]])
        abs.loc[len(abs)] = [t, avg, stdev, "Absorbance", "Community 2"]

        avg = np.average(data.loc["A", [7, 8, 9]])
        stdev = np.std(data.loc["A", [7, 8, 9]])
        abs.loc[len(abs)] = [t, avg, stdev, "Absorbance", "Community 3"]
    abs = abs.sort_values(by="exp_time")

    mcherry = pd.DataFrame(columns=["exp_time", "average", "stdev", "read", "species"])
    fs = listdir(join(dir, "at_oa/240921_ct_mono_batch", "plate_reader", "mCherry"))
    for f in sorted(fs):
        t = float(f[:-5])
        data = pd.read_excel(
            join(dir, "at_oa/240921_ct_mono_batch", "plate_reader", "mCherry", f),
            index_col=0,
        )
        avg = np.average(data.loc["A", [1, 2, 3]])
        stdev = np.std(data.loc["A", [1, 2, 3]])
        mcherry.loc[len(mcherry)] = [t, avg, stdev, "Fluorescence", "Community 1"]

        avg = np.average(data.loc["A", [4, 5, 6]])
        stdev = np.std(data.loc["A", [4, 5, 6]])
        mcherry.loc[len(mcherry)] = [t, avg, stdev, "Fluorescence", "Community 2"]

        avg = np.average(data.loc["A", [7, 8, 9]])
        stdev = np.std(data.loc["A", [7, 8, 9]])
        mcherry.loc[len(mcherry)] = [t, avg, stdev, "Fluorescence", "Community 3"]
    mcherry = mcherry.sort_values(by="exp_time")

    M0_abs = abs[abs["species"] == "Community 1"]
    M1_abs = abs[abs["species"] == "Community 2"]
    M2_abs = abs[abs["species"] == "Community 3"]

    M0_mcherry = mcherry[mcherry["species"] == "Community 1"]
    M1_mcherry = mcherry[mcherry["species"] == "Community 2"]
    M2_mcherry = mcherry[mcherry["species"] == "Community 3"]

    cfus = cfu_parser("at_oa/240921_ct_mono_batch")[0]
    ct_com1 = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "ct")]
    oa_com1 = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "oa")]
    ct_com2 = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "ct")]
    oa_com2 = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "oa")]
    ct_com3 = cfus[(cfus["reactor"] == "M2") & (cfus["species"] == "ct")]
    oa_com3 = cfus[(cfus["reactor"] == "M2") & (cfus["species"] == "oa")]

    return (
        df,
        M0,
        M1,
        M2,
        M0_abs,
        M1_abs,
        M2_abs,
        M0_mcherry,
        M1_mcherry,
        M2_mcherry,
        ct_com1,
        oa_com1,
        ct_com2,
        oa_com2,
        ct_com3,
        oa_com3,
    )


def plot_exp():
    (
        df,
        M0,
        M1,
        M2,
        M0_abs,
        M1_abs,
        M2_abs,
        M0_mcherry,
        M1_mcherry,
        M2_mcherry,
        at_com1,
        oa_com1,
        at_com2,
        oa_com2,
        at_com3,
        oa_com3,
    ) = parse_data()
    fig = make_subplots(
        rows=3,
        cols=3,
        shared_yaxes=False,
        shared_xaxes=True,
        row_titles=["OD", "mCherry", "mOrange"],
    )

    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["od_measured"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["od_measured"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["od_measured"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=1,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=2,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["FP3_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=3,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["FP3_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=3,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["FP3_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=3,
        col=3,
    )

    fig["layout"]["xaxis7"]["title"]["text"] = "Time [h]"
    fig["layout"]["xaxis8"]["title"]["text"] = "Time [h]"
    fig["layout"]["xaxis9"]["title"]["text"] = "Time [h]"
    yaxis = ["yaxis", "yaxis2", "yaxis3"]
    for yax in yaxis:
        fig["layout"][yax]["tick0"] = 0.1
        fig["layout"][yax]["dtick"] = 0.05
        fig["layout"][yax]["range"] = [0, 0.6]
        fig["layout"][yax]["tickmode"] = "linear"
    """yaxis = ["yaxis7", "yaxis8", "yaxis9"]
    for yax in yaxis:
        fig["layout"][yax]["type"] = "log"
"""
    return fig


def plot_triplicates():
    (
        df,
        M0,
        M1,
        M2,
        M0_abs,
        M1_abs,
        M2_abs,
        M0_mcherry,
        M1_mcherry,
        M2_mcherry,
        at_com1,
        oa_com1,
        at_com2,
        oa_com2,
        at_com3,
        oa_com3,
    ) = parse_data()
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_yaxes=True,
        shared_xaxes=True,
        row_titles=["OD", "mCherry"],
    )

    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["od_measured"],
            marker=dict(color="blue"),
            showlegend=True,
            name="Community 1",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["od_measured"],
            marker=dict(color="blue"),
            showlegend=True,
            line=dict(dash="dash"),
            name="Community 2",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["od_measured"],
            marker=dict(color="blue"),
            showlegend=True,
            line=dict(dash="dot"),
            name="Community 3",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
            line=dict(dash="dash"),
        ),
        row=2,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["raw_FP1_emit1"],
            marker=dict(color=colors["oa"]),
            showlegend=False,
            line=dict(dash="dot"),
        ),
        row=2,
        col=1,
    )

    fig = style_plot(fig, marker_size=3)
    fig.show()


fig = plot_exp()
fig.show()
