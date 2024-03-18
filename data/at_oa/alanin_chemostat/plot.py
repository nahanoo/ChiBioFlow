import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr
from plotly.subplots import make_subplots

df = fluorescence_paresr('at_oa/alanin_chemostat',
                         od_led_channel=['FP3_emit2'], od_led_converter=0.00001958)


fig = make_subplots(rows=3, cols=2, shared_yaxes=True, shared_xaxes=True, column_titles=[
                    'Community L-Alanine', 'Community L-Alanine + Thiamine'],
                    row_titles=['OD spectrometer', 'Density normalized GFP', 'Density normalized mCherry'])
ala = df[df['reactor'] == 'M0']
ala_t = df[df['reactor'] == 'M1']

fig.add_trace(go.Scatter(x=ala['exp_time'], y=ala['od_measured_led'], marker=dict(
    color='blue'),showlegend=False), row=1, col=1)
fig.add_trace(go.Scatter(x=ala_t['exp_time'], y=ala_t['od_measured_led'], marker=dict(
    color='blue'),showlegend=False), row=1, col=2)

fig.add_trace(go.Scatter(x=ala['exp_time'], y=ala['FP1_emit1'], marker=dict(
    color='#1B9E77'),showlegend=False), row=2, col=1)
fig.add_trace(go.Scatter(x=ala_t['exp_time'], y=ala_t['FP1_emit1'], marker=dict(
    color='#1B9E77'),showlegend=False), row=2, col=2)

fig.add_trace(go.Scatter(x=ala['exp_time'], y=ala['FP2_emit1'], marker=dict(
    color='#D95F02'),showlegend=False), row=3, col=1)
fig.add_trace(go.Scatter(x=ala_t['exp_time'], y=ala_t['FP2_emit1'], marker=dict(
    color='#D95F02'),showlegend=False), row=3, col=2)

fig['layout']['yaxis']['title'] = 'OD'
fig['layout']['yaxis3']['title'] = 'Fluorescence GFP'
fig['layout']['yaxis5']['title'] = 'Fluorescence mCherry'
fig['layout']['xaxis5']['title'] = 'Time [h]'
fig['layout']['xaxis6']['title'] = 'Time [h]'


