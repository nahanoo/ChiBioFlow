import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser, calibration
from plotly.subplots import make_subplots
from style_plot import *
import numpy as np
from models import *
from equations import *

params["v2_1"] = 0.23
params["K2_1"] = 0.02
params["q2_1"] = 0.329 / 7.5
params["v1_1"] = 0.5
params["K1_1"] = 0.2
params["q1_1"] = 0.281 / 7.5
params["M1"] = 7.5
params["D"] = 0.15
params["N01"] = 0.052 / 2
params["N02"] = 0.052 / 2

Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/241107_ct_oa")
Ms = [Ms[Ms["reactor"] == r] for r in ["M0", "M1", "M2"]]
# Ms = [M[M['exp_time'] >= 4] for M in Ms]
xs = np.linspace(0, Ms[0].iloc[-1]["exp_time"], int(Ms[0].iloc[-1]["exp_time"]) * 60)
xs = np.linspace(0, 40, 1000)
# y = species_12_niche_13(xs, f_dN1_1, f_dN2_13, f_dR12_1, f_dR2_3, params)
y = species_12_niche_1(xs, f_dN1_1, f_dN2_1, f_dR12_1, params)
k_eq, b_eq = calibration()
OD0 = 0.052
fig = make_subplots(
    rows=2,
    cols=1,
    shared_xaxes=True,
    subplot_titles=["OD", "Absolut abundances"],
    vertical_spacing=0.03,
)

Ms_restart = fluorescence_paresr("/home/eric/ChiBioFlow/data/241107_ct_oa_restart")
Ms_restart = [Ms_restart[Ms_restart["reactor"] == r] for r in ["M0", "M1", "M2"]]

for M, M_restart in zip(Ms, Ms_restart):
    M.index = range(len(M))
    # M.loc[:, 'exp_time'] = np.linspace(0, M.iloc[-1]['exp_time'], len(M))
    R0_raw = np.average(M.loc[:10]["od_measured"])
    t1 = M[(M["exp_time"] > 7.6) & (M["exp_time"] < 7.8)]
    R1_raw = np.average(t1["od_measured"])
    k_num = float(k_eq.subs({"OD1": OD0, "OD2": 0.3, R1: R0_raw, R2: R1_raw}))
    b_num = float(b_eq.subs({"OD1": OD0, "OD2": 0.3, R1: R0_raw, R2: R1_raw}))
    OD = k_num * np.log10(M["od_measured"]) + b_num
    fig.add_trace(
        go.Scatter(
            x=M["exp_time"][::10],
            y=OD[::10],
            marker=dict(color="#4F4F4F"),
            showlegend=False,
            hovertext=M["reactor"],
            name="OD",
        ),
        row=1,
        col=1,
    )
    """fig.add_trace(go.Scatter(x=M['exp_time'],
                             y=M['FP1_emit1'],
                             marker=dict(color=colors['oa']),
                             showlegend=False,
                             hovertext=M['reactor'],
                             name='OD'),
                  row=3,
                  col=1)"""
    t_end = M.iloc[-1]["exp_time"]
    M_restart.loc[:, "exp_time"] += 15
    OD = k_num * np.log10(M_restart["od_measured"]) + b_num
    fig.add_trace(
        go.Scatter(
            x=M_restart["exp_time"],
            y=OD,
            marker=dict(color="#4F4F4F"),
            showlegend=False,
            hovertext=M["reactor"],
            name="OD",
        ),
        row=1,
        col=1,
    )
"""    fig.add_trace(go.Scatter(x=M_restart['exp_time'],
                             y=M_restart['FP1_emit1'],
                             marker=dict(color=colors['oa']),
                             showlegend=False,
                             hovertext=M['reactor'],
                             name='OD'),
                  row=3,
                  col=1)"""
fig.add_trace(
    go.Scatter(
        x=xs,
        y=y[:, 0] + y[:, 1],
        line=dict(dash="dash"),
        marker=dict(color="#4F4F4F"),
        name="Model total",
    ),
    row=1,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=xs,
        y=y[:, 1],
        line=dict(dash="dash"),
        marker=dict(color=colors["oa"]),
        name="Model Oa",
    ),
    row=1,
    col=1,
)
fig.add_trace(
    go.Scatter(
        x=xs,
        y=y[:, 0],
        line=dict(dash="dash"),
        marker=dict(color=colors["ct"]),
        name="Model Ct",
    ),
    row=1,
    col=1,
)
"""fig.add_trace(go.Scatter(x=xs,
                         y=y[:, 1] / (y[:, 0] + y[:, 1]),
                         line=dict(dash='dash'),
                         marker=dict(color=colors['oa']),
                         showlegend=False,
                         name='Model Oa'),
              row=2,
              col=1)
fig.add_trace(go.Scatter(x=xs,
                         y=y[:, 0] / (y[:, 0] + y[:, 1]),
                         line=dict(dash='dash'),
                         marker=dict(color=colors['ct']),
                         showlegend=False,
                         name='Model Ct'),
              row=2,
              col=1)"""
cfus = cfu_parser("/home/eric/ChiBioFlow/data/at_oa/241107_ct_oa")[0]
Ms_cfus = [cfus[cfus["reactor"] == r] for r in sorted(set(cfus["reactor"]))]
Ms_ct = [M_cfu[M_cfu["species"] == "ct"] for M_cfu in Ms_cfus]
Ms_oa = [M_cfu[M_cfu["species"] == "oa"] for M_cfu in Ms_cfus]

for ct, oa in zip(Ms_ct, Ms_oa):
    fig.add_trace(
        go.Scatter(
            x=ct["sample_time"],
            y=ct["average"],
            line=dict(color=colors["ct"]),
            error_y=dict(type="data", array=ct["stdev"].to_list(), visible=True),
            name="Ct",
            showlegend=False,
            legend="legend2",
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=oa["sample_time"],
            y=oa["average"],
            line=dict(color=colors["oa"]),
            error_y=dict(type="data", array=oa["stdev"].to_list(), visible=True),
            name="Oa",
            showlegend=False,
            legend="legend2",
        ),
        row=2,
        col=1,
    )
"""    fig.add_trace(go.Scatter(x=ct['sample_time'],
                             y=ct['composition'],
                             line=dict(color=colors['ct']),
                             name='Ct',
                             showlegend=False,
                             legend='legend3'),
                  col=1,
                  row=3)"""
"""    fig.add_trace(go.Scatter(x=oa['sample_time'],
                             y=oa['composition'],
                             line=dict(color=colors['oa']),
                             name='Oa',
                             showlegend=False,
                             legend='legend3'),
                  col=1,
                  row=3)"""

fig["layout"]["yaxis1"]["title"] = "OD"
fig["layout"]["xaxis2"]["title"] = "Time [h]"
fig["layout"]["yaxis2"]["title"] = "CFUs/mL"

fig["data"][0]["showlegend"] = True
fig["data"][6]["showlegend"] = True
fig["data"][9]["showlegend"] = True
fig["data"][10]["showlegend"] = True

fig.update_layout(
    legend2={
        "y": 0.3,
    },
)
fig.update_layout(title="D = 0.15 1/h")
fig.update_layout(width=width / 2, height=height / 1.8)
fig = style_plot(fig, marker_size=1.5, left_margin=50)
fig["layout"]["yaxis2"]["type"] = "log"
# fig.show()

fig.write_image("plot.svg")
cfus = cfus[cfus["sample_time"] < 9]
cfus_ct = cfus[cfus["species"] == "ct"]
cfus_oa = cfus[cfus["species"] == "oa"]
ct0 = np.average(
    cfus[(cfus["sample_time"] == 0) & (cfus["species"] == "ct")]["average"]
)
ct0std = np.std(cfus[(cfus["sample_time"] == 0) & (cfus["species"] == "ct")]["average"])
ct1 = np.average(
    cfus[(cfus["sample_time"] == 7.8) & (cfus["species"] == "ct")]["average"]
)
ct1std = np.std(
    cfus[(cfus["sample_time"] == 7.8) & (cfus["species"] == "ct")]["average"]
)
oa0 = np.average(
    cfus[(cfus["sample_time"] == 0) & (cfus["species"] == "oa")]["average"]
)
oa0std = np.std(cfus[(cfus["sample_time"] == 0) & (cfus["species"] == "oa")]["average"])
oa1 = np.average(
    cfus[(cfus["sample_time"] == 7.8) & (cfus["species"] == "oa")]["average"]
)
oa1std = np.std(
    cfus[(cfus["sample_time"] == 7.8) & (cfus["species"] == "oa")]["average"]
)

dt = 7.8

vct = (np.log(ct1) - np.log(ct0)) / dt + 0.15
voa = (np.log(oa1) - np.log(oa0)) / dt + 0.15

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=[0, dt],
        y=[cfus_ct["average"].iloc[0], ct1],
        line=dict(color=colors["ct"]),
        error_y=dict(
            type="data", array=[cfus_ct["stdev"].iloc[0], ct1std], visible=True
        ),
        name="Ct",
    )
)
fig.add_trace(
    go.Scatter(
        x=[0, dt],
        y=[cfus_oa["average"].iloc[0], oa1],
        line=dict(color=colors["oa"]),
        error_y=dict(
            type="data", array=[cfus_oa["stdev"].iloc[0], oa1std], visible=True
        ),
        name="Oa",
    )
)
fig.update_layout(width=width / 2.5, height=height / 1.8 / 3)
fig["layout"]["yaxis"]["type"] = "log"
fig.update_xaxes(title="Time [h]")
fig.update_yaxes(title="CFUs/mL")
fig = style_plot(fig, marker_size=1.5, left_margin=60)

fig.write_image("zoomin.svg")


fig = go.Figure()
cfus = cfu_parser("/home/eric/ChiBioFlow/data/at_oa/241107_ct_oa")[0]
Ms_cfus = [cfus[cfus["reactor"] == r] for r in sorted(set(cfus["reactor"]))]
Ms_ct = [M_cfu[M_cfu["species"] == "ct"] for M_cfu in Ms_cfus]
Ms_oa = [M_cfu[M_cfu["species"] == "oa"] for M_cfu in Ms_cfus]

for ct, oa in zip(Ms_ct, Ms_oa):
    fig.add_trace(
        go.Scatter(
            x=ct["sample_time"],
            y=ct["average"],
            line=dict(color="#000080"),
            error_y=dict(type="data", array=ct["stdev"].to_list(), visible=True),
            name="#000080ff",
            showlegend=False,
            legend="legend2",
        ),
    )
    fig.add_trace(
        go.Scatter(
            x=oa["sample_time"],
            y=oa["average"],
            line=dict(color="#800080"),
            error_y=dict(type="data", array=oa["stdev"].to_list(), visible=True),
            name="#80008090",
            showlegend=False,
            legend="legend2",
        ),
    )
fig = style_plot(fig, font_size=10, marker_size=1.5)
fig.update_layout(width=200, height=170)
fig.update_layout(
    {
        "plot_bgcolor": "rgba(0, 0, 0, 0)",
        "paper_bgcolor": "rgba(0, 0, 0, 0)",
    },
)
fig.update_xaxes(title="Time", gridcolor="rgba(0, 0, 0, 0)")
fig.update_yaxes(title="CFUs/mL", type="log", gridcolor="rgba(0, 0, 0, 0)")
fig.write_image("plot_review.svg")
