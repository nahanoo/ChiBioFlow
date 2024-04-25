import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots
from math import log

df = fluorescence_paresr("at_oa/240425_oa_mono")
#cfus = cfu_parser("at_oa/at_mono")[0]

ala = df[df["reactor"] == "M0"]
ala_t = df[df["reactor"] == "M1"]
#at_ala = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "at")]
#at_ala_t = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "at")]


def plot_exp():
    fig = make_subplots(
        rows=2,
        cols=2,
        shared_yaxes=True,
        shared_xaxes=True,
        column_titles=["Oa Alanine", "Oa Alanine + Thiamine"],
        row_titles=["395 nm", "mCherry"],
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

    fig.show()


def max_growth_rate(D, t0, t1, N0,N1):
    return D + 1 / (t1 - t0) * log(N1 / N0)


D = 0.42
t0,t1 = 4.004, 6.09
N0,N1 = 5694,4928
u = max_growth_rate(D,t0,t1,N0,N1)