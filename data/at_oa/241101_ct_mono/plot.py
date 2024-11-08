import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, calibration
from plotly.subplots import make_subplots
from style_plot import *
import numpy as np
from models import *
from equations import *
from sympy import *

params['v1_1'] = 0.5
params['K1_1'] = 0.2
params['q1_1'] = 0.17 / 7.5
params['v2_1'] = 0.39
params['K2_1'] = 0.02
params['q2_1'] = 0.4 / 7.5
params['M1'] = 7.5
params['D'] = 0.28
params['N01'] = 0.05
params['N02'] = 0.04

Ms = fluorescence_paresr('/home/eric/ChiBioFlow/data/at_oa/241101_ct_mono')
Ms = [Ms[Ms['reactor'] == r] for r in ['M0', 'M1', 'M2']]
Ms = [M[M['exp_time'] < 23.5] for M in Ms]
xs = Ms[0]['exp_time']
y = species_1_niche_1(np.array(Ms[0]['exp_time']), f_dN1_1, f_dR1_1, params)
k_eq, b_eq = calibration()

fig = make_subplots(
    rows=1,
    cols=1,
    shared_xaxes=True,
    subplot_titles=['Total OD', 'Absolute abundances', 'Relative abundances'],
    vertical_spacing=0.03)
for M in Ms[1:]:
    M.index = range(len(M))
    R0_raw = np.average(M.loc[3:5]['od_measured'])
    R1_raw = np.average(M.loc[22.5 * 60:23.5 * 60]['od_measured'])
    k_num = float(
        k_eq.subs({
            'OD1': 0.05,
            'OD2': 0.141,
            R1: R0_raw,
            R2: R1_raw
        }))
    b_num = float(
        b_eq.subs({
            'OD1': 0.05,
            'OD2': 0.141,
            R1: R0_raw,
            R2: R1_raw
        }))
    OD = k_num * np.log10(M['od_measured']) + b_num
    fig.add_trace(go.Scatter(x=M['exp_time'],
                             y=OD,
                             marker=dict(color='#4F4F4F'),
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
fig["layout"]['yaxis1']["title"] = "OD"
fig["layout"]['xaxis1']["title"] = "Time [h]"
fig['data'][0]['showlegend'] = True
fig.update_layout(width=width, height=height)
fig = style_plot(fig)
fig.show()
fig.write_image('plot.png')
