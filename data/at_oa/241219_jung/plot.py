import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser, calibration
from plotly.subplots import make_subplots
from style_plot import *
import numpy as np
from equations import *

Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/241219_jung")
Ms = [Ms[Ms["reactor"] == r] for r in ["M0", "M1", "M2"]]
k_eq, b_eq = calibration()
OD0 = 0.16
fig = make_subplots(
    rows=4,
    cols=1,
    shared_xaxes=True,
    subplot_titles=["OD", "FP1", "FP2", "FP3"],
    vertical_spacing=0.03,
)
dfs = []
cols = [
    "exp_time",
    "od_measured",
    "raw_FP1_emit1",
    "raw_FP2_emit1",
    "raw_FP3_emit1",
    "reactor",
]
dfs = []

for M in Ms:
    M.index = range(len(M))
    # M.loc[:, 'exp_time'] = np.linspace(0, M.iloc[-1]['exp_time'], len(M))
    R0_raw = np.average(M.loc[:10]["od_measured"])
    t1 = M[(M["exp_time"] > 8.5) & (M["exp_time"] < 9)]
    R1_raw = np.average(t1["od_measured"])
    k_num = float(k_eq.subs({"OD1": OD0, "OD2": 3, R1: R0_raw, R2: R1_raw}))
    b_num = float(b_eq.subs({"OD1": OD0, "OD2": 3, R1: R0_raw, R2: R1_raw}))
    OD = k_num * np.log10(M["od_measured"]) + b_num
    M.loc[:, "od_measured"] = OD
    dfs.append(M[cols])

for M in Ms[:-1]:
    fig.add_trace(
        go.Scatter(
            x=M["exp_time"][::10],
            y=M["od_measured"][::10],
            marker=dict(color="#4F4F4F"),
            showlegend=False,
            hovertext=M["reactor"],
            name="OD",
        ),
        row=1,
        col=1,
    )
fig.add_trace(
    go.Scatter(
        x=Ms[0]["exp_time"][::10],
        y=Ms[0]["raw_FP1_emit1"][::10],
        marker=dict(color="#800080"),
        showlegend=False,
        hovertext=M["reactor"],
        name="FP1",
    ),
    row=2,
    col=1,
)

fig.add_trace(
    go.Scatter(
        x=Ms[0]["exp_time"][::10],
        y=Ms[0]["raw_FP2_emit1"][::10],
        marker=dict(color="#aa0000"),
        showlegend=False,
        hovertext=M["reactor"],
        name="FP1",
    ),
    row=3,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=Ms[0]["exp_time"][::10],
        y=Ms[0]["raw_FP3_emit1"][::10],
        marker=dict(color="#803300"),
        showlegend=False,
        hovertext=M["reactor"],
        name="FP1",
    ),
    row=4,
    col=1,
)

fig["layout"]["xaxis4"]["title"] = "Time [h]"
fig = style_plot(fig)
fig.update_layout(height=1000)
fig.write_image("plt.svg")

pd.concat(dfs).to_csv("data.csv", index=False)
