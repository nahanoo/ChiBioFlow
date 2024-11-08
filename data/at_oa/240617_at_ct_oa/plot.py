import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots
import pandas as pd
from os.path import join
from os import listdir, mkdir
import numpy as np
from glob import glob
from os.path import exists

colors = {
    "at": "#1B9E77",
    "oa": "#D95F02",
    "ct": "#7570B3",
    "<i>A. tumefaciens</i>": "#1B9E77",
    "<i>O. anthropi</i>": "#D95F02",
}
base_dir = join("/home", "eric", "ChiBioFlow", "data")

dirs = [join(base_dir, d) for d in ["240617_at_ct_oa", "240617_at_ct_oa_restart"]]
reactors = {"M0": [], "M1": [], "M2": [], "M3": []}

for dir in dirs:
    for reactor in reactors.keys():
        f = glob(join(dir, reactor, "*data.csv"))[0]
        reactors[reactor].append(pd.read_csv(f))

for reactor, dfs in reactors.items():
    trgt = join(base_dir, "at_oa", "240617_at_ct_oa", reactor)
    if not exists(trgt):
        mkdir(trgt)
    out = pd.concat(dfs)
    out["exp_time"] = np.linspace(0, len(out) * 60, len(out))
    out.to_csv(join(trgt, reactor + "_data.csv"), index=False)


df = fluorescence_paresr(
    join(base_dir, "at_oa/240617_at_ct_oa"),
    down_sample=False,
    window_size=20,
    od_filter=1,
)

M0 = df[df["reactor"] == "M0"]
M1 = df[df["reactor"] == "M1"]
M2 = df[df["reactor"] == "M2"]
M3 = df[df["reactor"] == "M3"]

fig = make_subplots(
    rows=3,
    cols=4,
    shared_yaxes=True,
    shared_xaxes=True,
    row_titles=["OD", "GFP", "mCherry"],
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
        x=M3["exp_time"],
        y=M3["od_measured"],
        marker=dict(color="blue"),
        showlegend=False,
    ),
    row=1,
    col=4,
)

fig.add_trace(
    go.Scatter(
        x=M0["exp_time"],
        y=M0["raw_FP1_emit1"],
        marker=dict(color=colors["ct"]),
        showlegend=False,
    ),
    row=2,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=M1["exp_time"],
        y=M1["raw_FP1_emit1"],
        marker=dict(color=colors["ct"]),
        showlegend=False,
    ),
    row=2,
    col=2,
)

fig.add_trace(
    go.Scatter(
        x=M2["exp_time"],
        y=M2["raw_FP1_emit1"],
        marker=dict(color=colors["at"]),
        showlegend=False,
    ),
    row=2,
    col=3,
)
fig.add_trace(
    go.Scatter(
        x=M3["exp_time"],
        y=M3["raw_FP1_emit1"],
        marker=dict(color=colors["at"]),
        showlegend=False,
    ),
    row=2,
    col=4,
)

fig.add_trace(
    go.Scatter(
        x=M0["exp_time"],
        y=M0["raw_FP2_emit1"],
        marker=dict(color=colors["oa"]),
        showlegend=False,
    ),
    row=3,
    col=1,
)

fig.add_trace(
    go.Scatter(
        x=M1["exp_time"],
        y=M1["raw_FP2_emit1"],
        marker=dict(color=colors["oa"]),
        showlegend=False,
    ),
    row=3,
    col=2,
)
fig.add_trace(
    go.Scatter(
        x=M1["exp_time"],
        y=M1["raw_FP2_emit1"],
        marker=dict(color=colors["oa"]),
        showlegend=False,
    ),
    row=3,
    col=3,
)
fig.add_trace(
    go.Scatter(
        x=M1["exp_time"],
        y=M1["raw_FP2_emit1"],
        marker=dict(color=colors["oa"]),
        showlegend=False,
    ),
    row=3,
    col=4,
)
