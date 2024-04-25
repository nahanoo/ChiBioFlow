import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots
from math import log

df = fluorescence_paresr("at_oa/240419_at")
cfus = cfu_parser("at_oa/240419_at")[0]

ala = df[df["reactor"] == "M0"]
ala_t = df[df["reactor"] == "M1"]
at_ala = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "at")]
at_ala_t = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "at")]


def plot_exp():
    fig = make_subplots(
        rows=3,
        cols=2,
        shared_yaxes=True,
        shared_xaxes=True,
        column_titles=["At Alanine", "At Alanine + Thiamine"],
        row_titles=["395 nm", "GFP", "CFUs/mL"],
    )

    fig.add_trace(
        go.Scatter(
            x=ala["exp_time"],
            y=ala["FP3_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=ala_t["exp_time"],
            y=ala_t["FP3_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=ala["exp_time"],
            y=ala["raw_FP1_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=ala_t["exp_time"],
            y=ala_t["raw_FP1_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=at_ala["sample_time"],
            y=at_ala["average"],
            marker=dict(color="#1B9E77"),
            error_y=dict(type="data", array=at_ala["stdev"].to_list(), visible=True),
            showlegend=False,
        ),
        row=3,
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
        row=3,
        col=2,
    )

    fig["layout"]["yaxis5"]["type"] = "log"
    fig["layout"]["yaxis6"]["type"] = "log"
    fig["layout"]["xaxis5"]["title"] = "Time [h]"
    fig["layout"]["xaxis6"]["title"] = "Time [h]"

    fig.show()


def max_growth_rate(D, t0, t1, N0,N1):
    return D + 1 / (t1 - t0) * log(N1 / N0)


def get_max_growth_rate():
    D = 0.26
    t0_i, t1_i = 5796, 6861
    t0_ala, t0_ala_t = df.loc[t0_i]["exp_time"].to_list()
    t1_ala, t1_ala_t = df.loc[t1_i]["exp_time"].to_list()

    N0_ala, N0_ala_t = df.loc[t0_i]["FP3_emit1"].to_list()
    N1_ala, N1_ala_t = df.loc[t1_i]["FP3_emit1"].to_list()

    u_ala = max_growth_rate(D, t0_ala, t1_ala, N0_ala, N1_ala)
    u_ala_t = max_growth_rate(D, t0_ala_t, t1_ala_t, N0_ala_t, N1_ala_t)

    return u_ala, u_ala_t

u_ala,u_ala_t = get_max_growth_rate()