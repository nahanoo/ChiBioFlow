import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser, calibration
from plotly.subplots import make_subplots
from style_plot import *
import numpy as np
from models import *
from equations import *

params['v2_1'] = 0.23
params['K2_1'] = 0.02
params['q2_1'] = 0.329 / 7.5
params['v1_1'] = 0.5
params['K1_1'] = 0.2
params['q1_1'] = 0.281 / 7.5
params['M1'] = 7.5
params['D'] = 0.15
params['N01'] = 0.052 / 2
params['N02'] = 0.052 / 2

Ms = fluorescence_paresr('/home/eric/ChiBioFlow/data/at_oa/241107_ct_oa')
Ms = [Ms[Ms['reactor'] == r] for r in ['M0', 'M1', 'M2']]
#Ms = [M[M['exp_time'] >= 4] for M in Ms]
xs = np.linspace(0, Ms[0].iloc[-1]['exp_time'],
                 int(Ms[0].iloc[-1]['exp_time']) * 60)
xs = np.linspace(0, 40, 1000)
#y = species_12_niche_13(xs, f_dN1_1, f_dN2_13, f_dR12_1, f_dR2_3, params)
y = species_12_niche_1(xs, f_dN1_1, f_dN2_1, f_dR12_1, params)
k_eq, b_eq = calibration()
OD0 = 0.052
fig = make_subplots(rows=3,
                    cols=1,
                    shared_xaxes=True,
                    subplot_titles=['Total OD', 'Relative abundances'],
                    vertical_spacing=0.03)

Ms_restart = fluorescence_paresr(
    '/home/eric/ChiBioFlow/data/241107_ct_oa_restart')
Ms_restart = [
    Ms_restart[Ms_restart['reactor'] == r] for r in ['M0', 'M1', 'M2']
]

for M, M_restart in zip(Ms, Ms_restart):
    M.index = range(len(M))
    #M.loc[:, 'exp_time'] = np.linspace(0, M.iloc[-1]['exp_time'], len(M))
    R0_raw = np.average(M.loc[:10]['od_measured'])
    t1 = M[(M['exp_time'] > 7.6) & (M['exp_time'] < 7.8)]
    R1_raw = np.average(t1['od_measured'])
    k_num = float(k_eq.subs({'OD1': OD0, 'OD2': 0.3, R1: R0_raw, R2: R1_raw}))
    b_num = float(b_eq.subs({'OD1': OD0, 'OD2': 0.3, R1: R0_raw, R2: R1_raw}))
    OD = k_num * np.log10(M['od_measured']) + b_num
    fig.add_trace(go.Scatter(x=M['exp_time'],
                             y=OD,
                             marker=dict(color='#4F4F4F'),
                             showlegend=False,
                             hovertext=M['reactor'],
                             name='OD'),
                  row=1,
                  col=1)
    fig.add_trace(go.Scatter(x=M['exp_time'],
                             y=M['FP1_emit1'],
                             marker=dict(color=colors['oa']),
                             showlegend=False,
                             hovertext=M['reactor'],
                             name='OD'),
                  row=3,
                  col=1)
    t_end = M.iloc[-1]['exp_time']
    M_restart.loc[:, 'exp_time'] += 15
    OD = k_num * np.log10(M_restart['od_measured']) + b_num
    fig.add_trace(go.Scatter(x=M_restart['exp_time'],
                             y=OD,
                             marker=dict(color='#4F4F4F'),
                             showlegend=False,
                             hovertext=M['reactor'],
                             name='OD'),
                  row=1,
                  col=1)
    fig.add_trace(go.Scatter(x=M_restart['exp_time'],
                             y=M_restart['raw_FP1_emit1'],
                             marker=dict(color=colors['oa']),
                             showlegend=False,
                             hovertext=M['reactor'],
                             name='OD'),
                  row=3,
                  col=1)
fig.add_trace(go.Scatter(x=xs,
                         y=y[:, 0] + y[:, 1],
                         line=dict(dash='dash'),
                         marker=dict(color='#4F4F4F'),
                         name='Model total'),
              row=1,
              col=1)
fig.add_trace(go.Scatter(x=xs,
                         y=y[:, 1],
                         line=dict(dash='dash'),
                         marker=dict(color=colors['oa']),
                         name='Model Oa'),
              row=1,
              col=1)
fig.add_trace(go.Scatter(x=xs,
                         y=y[:, 0],
                         line=dict(dash='dash'),
                         marker=dict(color=colors['ct']),
                         name='Model Ct'),
              row=1,
              col=1)
fig.add_trace(go.Scatter(x=xs,
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
              col=1)
fig["layout"]['yaxis1']["title"] = "OD"
fig["layout"]['xaxis1']["title"] = "Time [h]"
fig['data'][0]['showlegend'] = True
fig.update_layout(width=width, height=height * 1.5)
fig = style_plot(fig)
#fig.show()
#fig.write_image('plot.png')
