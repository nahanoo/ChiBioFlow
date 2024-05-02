import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots
from math import log
from sympy import symbols, solve, Eq
from sympy import init_printing

init_printing()

def parse_data():
    df = fluorescence_paresr("at_oa/240429_at_oa_mono")
    # cfus = cfu_parser("at_oa/at_mono")[0]
    df.index = range(len(df))
    to_pop = [
        122,
        1427,
        2991,
        1291,
        1292,
        1293,
        1294,
        1295,
        1296,
        1297,
        1298,
        4330,
        4331,
        4332,
        4333,
        4334,
        4335,
        4336,
    ]
    for i, OD in enumerate(df["od_measured"]):
        if OD < 0.2:
            to_pop.append(i)

    for i in to_pop:
        df = df.drop(i)
    M0 = df[df["reactor"] == "M0"]
    M0.loc[1299:, "od_measured"] = M0.loc[1299:, "od_measured"] - 0.19814100000000012
    M1 = df[df["reactor"] == "M1"]
    M1.loc[:4336, "od_measured"] = M1.loc[:4336, "od_measured"] + 0.28931299999999993 + 0.029867000000000088
    M2 = df[df["reactor"] == "M2"]
    return df,M0,M1,M2

# at_ala = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "at")]
# at_ala_t = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "at")]


def plot_exp():
    df,M0,M1,M2 = parse_data()
    fig = make_subplots(
        rows=3,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        column_titles=["At Alanine", "At Alanine + Thiamine, Oa Alanine + Thiamine"],
        row_titles=["OD", "Fluorescence",'Density LED 395 nm'],
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
            x=M1["exp_time"],
            y=M1["FP2_emit1"],
            marker=dict(color="blue"),
            showlegend=False,
        ),
        row=3,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=M2["exp_time"],
            y=M2["FP2_emit1"],
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
    fig.update_xaxes(title='Time [h]')

    fig.show()


#plot_exp()


