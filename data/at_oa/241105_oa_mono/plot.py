import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser, calibration
from plotly.subplots import make_subplots
from style_plot import *
import numpy as np
from models import *
from equations import *


def plot_oa_chibio():
    params['v1_1'] = 0.23
    params['K1_1'] = 0.02
    params['q1_1'] = 0.329 / 7.5
    params['v2_1'] = 0.39
    params['K2_1'] = 0.02
    params['q2_1'] = 0.4 / 7.5
    params['M1'] = 7.5
    params['D'] = 0.14
    params['N01'] = 0.027
    params['N02'] = 0.04

    Ms = fluorescence_paresr('/home/eric/ChiBioFlow/data/at_oa/241105_oa_mono')
    Ms = [Ms[Ms['reactor'] == r] for r in ['M0', 'M2']]
    Ms = [M[M['exp_time'] >= 4] for M in Ms]
    xs = np.linspace(0, Ms[0].iloc[-1]['exp_time'],
                     int(Ms[0].iloc[-1]['exp_time']) * 60)
    y = species_1_niche_1(xs, f_dN1_1, f_dR1_1, params)
    k_eq, b_eq = calibration()
    OD0 = 0.027
    fig = make_subplots(rows=1,
                        cols=1,
                        shared_xaxes=True,
                        vertical_spacing=0.03)
    for M in Ms:
        M.index = range(len(M))
        M.loc[:, 'exp_time'] = np.linspace(0, M.iloc[-1]['exp_time'], len(M))
        t_R1 = M[(M['exp_time'] < 23.5) & (M['exp_time'] < 24)]
        R0_raw = np.average(M.loc[:10]['od_measured'])
        R1_raw = np.average(M.iloc[-60:]['od_measured'])
        #R1_raw = np.average(t_R1['od_measured'])
        k_num = float(
            k_eq.subs({
                'OD1': OD0,
                'OD2': 0.3,
                R1: R0_raw,
                R2: R1_raw
            }))
        b_num = float(
            b_eq.subs({
                'OD1': OD0,
                'OD2': 0.3,
                R1: R0_raw,
                R2: R1_raw
            }))
        OD = k_num * np.log10(M['od_measured']) + b_num
        fig.add_trace(go.Scatter(x=M['exp_time'][::10],
                                 y=OD[::10],
                                 marker=dict(color=colors['oa']),
                                 showlegend=False,
                                 hovertext=M['reactor'],
                                 name='OD'),
                      row=1,
                      col=1)
    fig.add_trace(go.Scatter(x=xs,
                             y=y[:, 0],
                             line=dict(dash='dash'),
                             marker=dict(color='#4F4F4F'),
                             name='Model'),
                  row=1,
                  col=1)
    b = 0.3 / np.average(Ms[0].iloc[-5:]['FP3_emit1'])
    #fig.add_trace(go.Scatter(x=Ms[0]['exp_time'], y=Ms[0]['FP3_emit1'] * b))
    ##b = OD0 / np.average(Ms[2].loc[:3]['FP3_emit1'])
    #fig.add_trace(
    #    go.Scatter(x=Ms[0]['exp_time'],
    #               y=[get_od(raw) for raw in Ms[0]['od_measured']]))
    fig["layout"]['yaxis1']["title"] = "OD"
    fig["layout"]['xaxis1']["title"] = "Time [h]"
    fig['data'][0]['showlegend'] = True
    fig.update_layout(width=width, height=height)
    fig = style_plot(fig)
    #fig.show()
    fig.write_image('plot.png')
    return fig


fig = plot_oa_chibio()
fig.show()
