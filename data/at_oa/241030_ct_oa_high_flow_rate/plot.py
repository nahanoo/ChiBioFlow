import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser
from plotly.subplots import make_subplots
from style_plot import *
import numpy as np
from models import *
from equations import *

params['v1_1'] = 0.39
params['K1_1'] = 0.2
params['q1_1'] = 0.4 / 7.5
params['v2_1'] = 0.39
params['K2_1'] = 0.02
params['q2_1'] = 0.4 / 7.5
params['M1'] = 7.5
params['D'] = 0
params['N01'] = 0.04
params['N02'] = 0.04

xs = np.linspace(0, 60, 1000)
y = species_12_niche_13(xs, f_dN1_1, f_dN2_13, f_dR12_1, f_dR2_3, params)

Ms, cfus = ct_oa_parser(
    '/home/eric/ChiBioFlow/data/at_oa/241030_ct_oa_high_flow_rate')

fig = make_subplots(
    rows=2,
    cols=1,
    shared_xaxes=True,
    subplot_titles=['Total OD', 'Absolute abundances', 'Relative abundances'],
    vertical_spacing=0.03)
for M in Ms:
    M = M[M['exp_time'] < 30]
    b = 0.043 / M['FP3_emit1'][0]
    OD_FP3 = M['FP3_emit1'] * b
    if M.iloc[0]['reactor'] == 'M0':
        y = OD_FP3
    else:
        y = M['od_measured']

    fig.add_trace(go.Scatter(x=M['exp_time'],
                             y=y,
                             marker=dict(color='#4F4F4F'),
                             showlegend=False,
                             hovertext=M['reactor']),
                  row=1,
                  col=1)
"""fig.add_trace(
    go.Scatter(x=xs,
               y=y[:, 0] + y[:, 1],
               marker=dict(color='#4F4F4F'),
               line=dict(dash='dash')))
fig.add_trace(
    go.Scatter(x=xs,
               y=y[:, 0],
               marker=dict(color=colors['ct']),
               line=dict(dash='dash')))
fig.add_trace(
    go.Scatter(x=xs,
               y=y[:, 1],
               marker=dict(color=colors['oa']),
               line=dict(dash='dash')))"""
for species in ['ct', 'oa']:
    for M in cfus[species]:
        M.loc[:, 'average'] = M['average'].replace(0, np.nan)
        M = M[M['sample_time'] < 30]
        fig.add_trace(go.Scatter(x=M['sample_time'],
                                 y=M['average'],
                                 name=abb[species],
                                 marker=dict(color=colors[species]),
                                 error_y=dict(type="data",
                                              array=M["stdev"].to_list(),
                                              visible=True),
                                 showlegend=False,
                                 hovertext=M['reactor']),
                      row=2,
                      col=1)
"""
for species in ['ct', 'oa']:
    for M in cfus[species]:
        fig.add_trace(go.Scatter(x=M['sample_time'],
                                 y=M['composition'],
                                 marker=dict(color=colors[species]),
                                 showlegend=False,
                                 hovertext=M['reactor']),
                      row=3,
                      col=1)"""

fig['data'][3]['showlegend'] = True
fig['data'][6]['showlegend'] = True
fig["layout"]['yaxis2']["type"] = "log"
fig["layout"]['yaxis1']["title"] = "OD"
fig["layout"]['yaxis2']["title"] = "CFUs/mL"
fig["layout"]['xaxis2']["title"] = "Time [h]"
fig.update_layout(
    legend=dict(
        title='Species',
        y=0.6,  # Position the legend vertically; adjust as needed for placement
    ),
    width=width,
    height=height)
fig = style_plot(fig)
fig.show()
