import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots

df = fluorescence_paresr(
    "at_oa/240315_at_oa", od_led_channel=["FP3_emit2"], od_led_converter=0.00001958
)

cfus = cfu_parser("at_oa/240315_at_oa")[0]
fig = make_subplots(
    rows=4,
    cols=2,
    shared_yaxes=True,
    shared_xaxes=True,
    column_titles=["Community L-Alanine", "Community L-Alanine + Thiamine"],
    row_titles=[
        "OD spectrometer",
        "Density normalized GFP",
        "Density normalized mCherry",
    ],
)

ala = df[df["reactor"] == "M0"]
ala_t = df[df["reactor"] == "M1"]

at_ala = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "at")]
at_ala_t = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "at")]

oa_ala = cfus[(cfus["reactor"] == "M0") & (cfus["species"] == "oa")]
oa_ala_t = cfus[(cfus["reactor"] == "M1") & (cfus["species"] == "oa")]


fig.add_trace(
    go.Scatter(
        x=ala["exp_time"],
        y=ala["od_measured"],
        marker=dict(color="blue"),
        showlegend=False,
    ),
    row=1,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=ala_t["exp_time"],
        y=ala_t["od_measured"],
        marker=dict(color="blue"),
        showlegend=False,
    ),
    row=1,
    col=2,
)

fig.add_trace(
    go.Scatter(
        x=ala["exp_time"],
        y=ala["FP1_emit1"],
        marker=dict(color="#1B9E77"),
        showlegend=True,
        name="At",
    ),
    row=2,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=ala_t["exp_time"],
        y=ala_t["FP1_emit1"],
        marker=dict(color="#1B9E77"),
        showlegend=False,
    ),
    row=2,
    col=2,
)

fig.add_trace(
    go.Scatter(
        x=ala["exp_time"],
        y=ala["FP2_emit1"],
        marker=dict(color="#D95F02"),
        showlegend=True,
        name="Oa",
    ),
    row=2,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=ala_t["exp_time"],
        y=ala_t["FP2_emit1"],
        marker=dict(color="#D95F02"),
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

fig.add_trace(
    go.Scatter(
        x=oa_ala["sample_time"],
        y=oa_ala["average"],
        marker=dict(color="#D95F02"),
        error_y=dict(type="data", array=oa_ala["stdev"].to_list(), visible=True),
        showlegend=False,
    ),
    row=3,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=oa_ala_t["sample_time"],
        y=oa_ala_t["average"],
        marker=dict(color="#D95F02"),
        error_y=dict(type="data", array=oa_ala_t["stdev"].to_list(), visible=True),
        showlegend=False,
    ),
    row=3,
    col=2,
)

fig.add_trace(
    go.Scatter(
        x=at_ala["sample_time"],
        y=at_ala["composition"],
        marker=dict(color="#1B9E77"),
        showlegend=False,
    ),
    row=4,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=at_ala_t["sample_time"],
        y=at_ala_t["composition"],
        marker=dict(color="#1B9E77"),
        showlegend=False,
    ),
    row=4,
    col=2,
)

fig.add_trace(
    go.Scatter(
        x=oa_ala["sample_time"],
        y=oa_ala["composition"],
        marker=dict(color="#D95F02"),
        showlegend=False,
    ),
    row=4,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=oa_ala_t["sample_time"],
        y=oa_ala_t["composition"],
        marker=dict(color="#D95F02"),
        showlegend=False,
    ),
    row=4,
    col=2,
)


fig["layout"]["yaxis"]["title"] = "OD"
fig["layout"]["yaxis3"]["title"] = "Fluorescence GFP"
fig["layout"]["yaxis5"]["title"] = "CFUs/mL"
fig["layout"]["yaxis7"]["title"] = "Composition"
fig["layout"]["xaxis7"]["title"] = "Time [h]"
fig["layout"]["xaxis8"]["title"] = "Time [h]"
fig["layout"]["yaxis5"]["type"] = "log"
fig["layout"]["yaxis6"]["type"] = "log"
