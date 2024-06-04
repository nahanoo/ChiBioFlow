import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots
from math import log
from sympy import symbols, solve, Eq
from sympy import init_printing
import pandas as pd
from os.path import join
from os import listdir
import numpy as np

init_printing()


def parse_data():
    dir = join("/", "home", "eric", "ChiBioFlow", "data")
    df = fluorescence_paresr(join(dir, "at_oa/240429_at_oa_mono"))

    restart = fluorescence_paresr(join(dir, "240506_at_oa_mono_restart"))
    restart["exp_time"] += 77.16144722222224
    df = pd.concat([df, restart])
    df.index = range(len(df))
    to_pop = []
    for i, OD, time in zip(df.index, df["od_measured"], df["exp_time"]):
        if (OD < 0.2) & (time < 85):
            to_pop.append(i)
    for i in to_pop:
        df = df.drop(i)
    M0 = df[df["reactor"] == "M0"]
    M1 = df[df["reactor"] == "M1"]
    M2 = df[df["reactor"] == "M2"]

    abs = pd.DataFrame(columns=["exp_time", "average", "stdev", "read", "species"])
    fs = listdir(join(dir, "at_oa/240429_at_oa_mono", "plate_reader", "Absorbance"))
    for f in sorted(fs):
        t = float(f[:-5])
        data = pd.read_excel(
            join(dir, "at_oa/240429_at_oa_mono", "plate_reader", "Absorbance", f),
            index_col=0,
        )
        avg = np.average(data.loc["A", [1, 2, 3]])
        stdev = np.std(data.loc["A", [1, 2, 3]])
        abs.loc[len(abs)] = [t, avg, stdev, "Absorbance", "At Alanine"]

        avg = np.average(data.loc["A", [4, 5, 6]])
        stdev = np.std(data.loc["A", [4, 5, 6]])
        abs.loc[len(abs)] = [t, avg, stdev, "Absorbance", "At Alanine + Thiamine"]

        avg = np.average(data.loc["A", [7, 8, 9]])
        stdev = np.std(data.loc["A", [7, 8, 9]])
        abs.loc[len(abs)] = [t, avg, stdev, "Absorbance", "Oa Alanine + Thiamine"]
    abs = abs.sort_values(by="exp_time")

    gfp = pd.DataFrame(columns=["exp_time", "average", "stdev", "read", "species"])
    fs = listdir(join(dir, "at_oa/240429_at_oa_mono", "plate_reader", "GFP"))
    for f in sorted(fs):
        t = float(f[:-5])
        data = pd.read_excel(
            join(dir, "at_oa/240429_at_oa_mono", "plate_reader", "GFP", f), index_col=0
        )
        avg = np.average(data.loc["A", [1, 2, 3]])
        stdev = np.std(data.loc["A", [1, 2, 3]])
        gfp.loc[len(gfp)] = [t, avg, stdev, "Absorbance", "At Alanine"]

        avg = np.average(data.loc["A", [4, 5, 6]])
        stdev = np.std(data.loc["A", [4, 5, 6]])
        gfp.loc[len(gfp)] = [t, avg, stdev, "Absorbance", "At Alanine + Thiamine"]
    gfp = gfp.sort_values(by="exp_time")

    mcherry = pd.DataFrame(columns=["exp_time", "average", "stdev", "read", "species"])
    fs = listdir(join(dir, "at_oa/240429_at_oa_mono", "plate_reader", "mCherry"))
    for f in sorted(fs):
        t = float(f[:-5])
        data = pd.read_excel(
            join(dir, "at_oa/240429_at_oa_mono", "plate_reader", "mCherry", f),
            index_col=0,
        )
        avg = np.average(data.loc["A", [7, 8, 9]])
        stdev = np.std(data.loc["A", [7, 8, 9]])
        mcherry.loc[len(mcherry)] = [
            t,
            avg,
            stdev,
            "Absorbance",
            "Oa Alanine + Thiamine",
        ]
    mcherry = mcherry.sort_values(by="exp_time")

    M0_abs = abs[abs["species"] == "At Alanine"]
    M1_abs = abs[abs["species"] == "At Alanine + Thiamine"]
    M2_abs = abs[abs["species"] == "Oa Alanine + Thiamine"]

    M0_gfp = gfp[gfp["species"] == "At Alanine"]
    M1_gfp = gfp[gfp["species"] == "At Alanine + Thiamine"]
    M2_mcherry = mcherry[mcherry["species"] == "Oa Alanine + Thiamine"]

    cfus = cfu_parser("at_oa/240429_at_oa_mono")[0]
    at_ala = cfus[cfus["reactor"] == "M0"]
    at_ala_t = cfus[cfus["reactor"] == "M1"]
    oa_ala_t = cfus[cfus["reactor"] == "M2"]

    return (
        df,
        M0,
        M1,
        M2,
        M0_abs,
        M1_abs,
        M2_abs,
        M0_gfp,
        M1_gfp,
        M2_mcherry,
        at_ala,
        at_ala_t,
        oa_ala_t,
    )


# at_ala = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "at")]
# at_ala_t = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "at")]


def plot_exp():
    (
        df,
        M0,
        M1,
        M2,
        M0_abs,
        M1_abs,
        M2_abs,
        M0_gfp,
        M1_gfp,
        M2_mcherry,
        at_ala,
        at_ala_t,
        oa_ala_t,
    ) = parse_data()
    fig = make_subplots(
        rows=3,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        column_titles=["At Alanine", "At Alanine + Thiamine", " Oa Alanine + Thiamine"],
        row_titles=["OD", "Fluorescence", "Density LED 395 nm"],
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
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["raw_FP1_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["raw_FP1_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=2,
        col=3,
    )
    fig.add_trace(
        go.Scatter(
            x=M0["exp_time"],
            y=M0["FP2_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=3,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=M1["exp_time"],
            y=M1["FP2_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=3,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["FP2_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=3,
        col=3,
    )
    fig.update_xaxes(title="Time [h]")

    fig.show()
    # fig.write_html('mono_culttures.html')


def plot_plate_reader():

    (
        df,
        M0,
        M1,
        M2,
        M0_abs,
        M1_abs,
        M2_abs,
        M0_gfp,
        M1_gfp,
        M2_mcherry,
        at_ala,
        at_ala_t,
        oa_ala_t,
    ) = parse_data()
    fig = make_subplots(
        rows=4,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        column_titles=["At Alanine", "At Alanine + Thiamine", "Oa Alanine + Thiamine"],
        row_titles=["OD Chi.Bio", "OD plate reader", "Fluorescence plate reader"],
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
            x=M0_abs["exp_time"],
            y=M0_abs["average"],
            marker=dict(color="blue"),
            showlegend=False,
            error_y=dict(type="data", array=M0_abs["stdev"].to_list(), visible=True),
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1_abs["exp_time"],
            y=M1_abs["average"],
            marker=dict(color="blue"),
            showlegend=False,
            error_y=dict(type="data", array=M1_abs["stdev"].to_list(), visible=True),
        ),
        row=2,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=M2_abs["exp_time"],
            y=M2_abs["average"],
            marker=dict(color="blue"),
            showlegend=False,
            error_y=dict(type="data", array=M2_abs["stdev"].to_list(), visible=True),
        ),
        row=2,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=M0_gfp["exp_time"],
            y=M0_gfp["average"],
            marker=dict(color="blue"),
            showlegend=False,
            error_y=dict(type="data", array=M0_gfp["stdev"].to_list(), visible=True),
        ),
        row=3,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=M1_gfp["exp_time"],
            y=M1_gfp["average"],
            marker=dict(color="blue"),
            showlegend=False,
            error_y=dict(type="data", array=M1_gfp["stdev"].to_list(), visible=True),
        ),
        row=3,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=M2_mcherry["exp_time"],
            y=M2_mcherry["average"],
            marker=dict(color="blue"),
            showlegend=False,
            error_y=dict(
                type="data", array=M2_mcherry["stdev"].to_list(), visible=True
            ),
        ),
        row=3,
        col=3,
    )

    fig.add_trace(
        go.Scatter(
            x=at_ala["sample_time"],
            y=at_ala["average"],
            marker=dict(color="#1B9E77"),
            error_y=dict(type="data", array=at_ala["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=4,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=at_ala_t["sample_time"],
            y=at_ala_t["average"],
            marker=dict(color="#1B9E77"),
            error_y=dict(type="data", array=at_ala_t["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=4,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=oa_ala_t["sample_time"],
            y=oa_ala_t["average"],
            marker=dict(color="#D95F02"),
            error_y=dict(type="data", array=oa_ala_t["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=4,
        col=3,
    )

    # fig.for_each_xaxis(lambda x: x.update(showgrid=True, zeroline=False, range=[0, 60]))
    fig["layout"]["yaxis10"]["type"] = "log"
    fig["layout"]["yaxis11"]["type"] = "log"
    fig["layout"]["yaxis12"]["type"] = "log"

    fig.show()


