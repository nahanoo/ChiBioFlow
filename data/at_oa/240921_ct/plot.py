import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser
from plotly.subplots import make_subplots
from style_plot import *
from models import *
from equations import *
import numpy as np

params['D'] = 0.1
params['K2_1'] = 0.45
params['v2_1'] = 0.45
params['q2_1'] = 0.25 / 7.5
params['M1'] = 7.5
#params['N02'] = 0.1

Ms = fluorescence_paresr('/home/eric/ChiBioFlow/data/at_oa/240921_ct')
Ms = [Ms[Ms['reactor'] == i] for i in ['M0', 'M1', 'M2']]
M1 = Ms[1]
M1.index = range(len(M1))
for i, t in enumerate(M1['exp_time']):
    if t > 6.6:
        break
M1.loc[i + 5:, 'od_measured'] -= 0.22
Ms[1] = M1
for i, M in enumerate(Ms):
    M = M[(M['exp_time'] > 1) & (M['exp_time'] < 10)]
    raw = M[M['exp_time'] > 1].iloc[0]['FP3_emit1']
    b = raw / 0.05
    M.loc[:, 'FP3_emit1'] = M['FP3_emit1'] / b
    Ms[i] = M

fig = make_subplots(rows=1, cols=1, shared_xaxes=True, subplot_titles=['OD'])

xs = np.linspace(0, 10)
y = species_2_niche_1(xs, f_dN2_1, f_dR2_1, params)
legend = True
for M in Ms:
    fig.add_trace(go.Scatter(x=M['exp_time'][::10],
                             y=M['od_measured'][::10],
                             marker=dict(color='#4F4F4F'),
                             showlegend=legend,
                             name='OD',
                             hovertext=M['reactor']),
                  row=1,
                  col=1)
    legend = False
fig.add_trace(
    go.Scatter(x=xs,
               y=y[:, 0],
               line=dict(color='#4F4F4F', dash='dash'),
               name='Model'))

fig["layout"]['yaxis1']["title"] = "OD"
fig["layout"]['xaxis1']["title"] = "Time [h]"
fig.update_layout(width=width / 2.8,
                  height=height / 3.5,
                  title='D = 0.1 [1/h]')
fig = style_plot(fig, marker_size=1.5)
fig.write_image('plot.svg')
