import pandas as pd
import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser
from plotly.subplots import make_subplots
from os.path import join


def plot_exp():
    df = fluorescence_paresr('at_oa/240412_at_oa')
    cfus = cfu_parser('at_oa/240412_at_oa')[0]

    com_ala = df[df['reactor'] == 'M0']
    com_ala_t = df[df['reactor'] == 'M1']

    at_ala = cfus[(cfus['reactor'] == 'M0') & (cfus['species'] == 'at')]
    at_ala_t = cfus[(cfus['reactor'] == 'M1') & (cfus['species'] == 'at')]

    oa_ala = cfus[(cfus['reactor'] == 'M0') & (cfus['species'] == 'oa')]
    oa_ala_t = cfus[(cfus['reactor'] == 'M1') & (cfus['species'] == 'oa')]

    fig = make_subplots(rows=5, cols=2, shared_yaxes=False, shared_xaxes=True, column_titles=[
                        'Com. Alanine', 'Com. Alanine + Thiamine'],
                        row_titles=['OD', 'GFP', 'mCherry', 'CFUs/mL', 'Composition'])

    fig.add_trace(go.Scatter(x=com_ala['exp_time'], y=com_ala['od_measured'], marker=dict(
        color='blue'), showlegend=False), row=1, col=1)
    fig.add_trace(go.Scatter(x=com_ala['exp_time'], y=com_ala['raw_FP1_emit1'], marker=dict(
        color='blue'), showlegend=False), row=2, col=1)
    fig.add_trace(go.Scatter(x=com_ala['exp_time'], y=com_ala['raw_FP2_emit1'], marker=dict(
        color='blue'), showlegend=False), row=3, col=1)

    fig.add_trace(go.Scatter(x=com_ala_t['exp_time'], y=com_ala_t['od_measured'], marker=dict(
        color='blue'), showlegend=False), row=1, col=2)
    fig.add_trace(go.Scatter(x=com_ala_t['exp_time'], y=com_ala_t['raw_FP1_emit1'], marker=dict(
        color='blue'), showlegend=False), row=2, col=2)
    fig.add_trace(go.Scatter(x=com_ala_t['exp_time'], y=com_ala_t['raw_FP2_emit1'], marker=dict(
        color='blue'), showlegend=False), row=3, col=2)

    fig.add_trace(go.Scatter(x=at_ala['sample_time'], y=at_ala['average'],
                             marker=dict(color='#1B9E77'), error_y=dict(type='data', array=at_ala['stdev'].to_list(), visible=True), showlegend=False), row=4, col=1)
    fig.add_trace(go.Scatter(x=at_ala_t['sample_time'], y=at_ala_t['average'],
                             marker=dict(color='#1B9E77'), error_y=dict(type='data', array=at_ala_t['stdev'].to_list(), visible=True), showlegend=False), row=4, col=2)

    fig.add_trace(go.Scatter(x=oa_ala['sample_time'], y=oa_ala['average'],
                             marker=dict(color='#D95F02'), error_y=dict(type='data', array=oa_ala['stdev'].to_list(), visible=True), showlegend=False), row=4, col=1)
    fig.add_trace(go.Scatter(x=oa_ala_t['sample_time'], y=oa_ala_t['average'],
                             marker=dict(color='#D95F02'), error_y=dict(type='data', array=oa_ala_t['stdev'].to_list(), visible=True), showlegend=False), row=4, col=2)

    fig.add_trace(go.Scatter(x=at_ala['sample_time'], y=at_ala['composition'],
                             marker=dict(color='#1B9E77'), showlegend=False), row=5, col=1)
    fig.add_trace(go.Scatter(x=at_ala_t['sample_time'], y=at_ala_t['composition'],
                             marker=dict(color='#1B9E77'),  showlegend=False), row=5, col=2)

    fig.add_trace(go.Scatter(x=oa_ala['sample_time'], y=oa_ala['composition'],
                             marker=dict(color='#D95F02'), showlegend=False), row=5, col=1)
    fig.add_trace(go.Scatter(x=oa_ala_t['sample_time'], y=oa_ala_t['composition'],
                             marker=dict(color='#D95F02'), showlegend=False), row=5, col=2)

    fig['layout']['yaxis7']['type'] = 'log'
    fig['layout']['yaxis8']['type'] = 'log'
    fig['layout']['xaxis9']['title'] = 'Time [h]'
    fig['layout']['xaxis10']['title'] = 'Time [h]'

    return fig

fig = plot_exp()
fig.show()

def plate_reader():
    ala_plate = pd.read_csv(join('M0', 'plate_reader.csv'))
    ala_t_plate = pd.read_csv(join('M1', 'plate_reader.csv'))

    df = fluorescence_paresr('at_oa/240412_at_oa')


    com_ala = df[df['reactor'] == 'M0']
    com_ala_t = df[df['reactor'] == 'M1']

    cfus = cfu_parser('at_oa/240412_at_oa')[0]

    at_ala = cfus[(cfus['reactor'] == 'M0') & (cfus['species'] == 'at')]
    at_ala_t = cfus[(cfus['reactor'] == 'M1') & (cfus['species'] == 'at')]

    oa_ala = cfus[(cfus['reactor'] == 'M0') & (cfus['species'] == 'oa')]
    oa_ala_t = cfus[(cfus['reactor'] == 'M1') & (cfus['species'] == 'oa')]

    fig = make_subplots(rows=7, cols=2, shared_yaxes=False, shared_xaxes=True, column_titles=[
                        'Com. Alanine', 'Com. Alanine + Thiamine'],
                        row_titles=['OD ChiBio', 'GFP ChiBio', 'mCherry ChiBio','OD PR','GFP PR','mCherry PR', 'CFUs/mL'])


    fig.add_trace(go.Scatter(x=com_ala['exp_time'], y=com_ala['od_measured'], marker=dict(
        color='blue'), showlegend=False), row=1, col=1)
    fig.add_trace(go.Scatter(x=com_ala['exp_time'], y=com_ala['raw_FP1_emit1'], marker=dict(
        color='blue'), showlegend=False), row=2, col=1)
    fig.add_trace(go.Scatter(x=com_ala['exp_time'], y=com_ala['raw_FP2_emit1'], marker=dict(
        color='blue'), showlegend=False), row=3, col=1)

    fig.add_trace(go.Scatter(x=com_ala_t['exp_time'], y=com_ala_t['od_measured'], marker=dict(
        color='blue'), showlegend=False), row=1, col=2)
    fig.add_trace(go.Scatter(x=com_ala_t['exp_time'], y=com_ala_t['raw_FP1_emit1'], marker=dict(
        color='blue'), showlegend=False), row=2, col=2)
    fig.add_trace(go.Scatter(x=com_ala_t['exp_time'], y=com_ala_t['raw_FP2_emit1'], marker=dict(
        color='blue'), showlegend=False), row=3, col=2)

    fig.add_trace(go.Scatter(x=ala_plate['time'], y=ala_plate['OD'], marker=dict(
        color='blue'), showlegend=False), row=4, col=1)
    fig.add_trace(go.Scatter(x=ala_t_plate['time'], y=ala_t_plate['OD'], marker=dict(
        color='blue'), showlegend=False), row=4, col=2)

    fig.add_trace(go.Scatter(x=ala_plate['time'], y=ala_plate['GFP'], marker=dict(
        color='blue'), showlegend=False), row=5, col=1)
    fig.add_trace(go.Scatter(x=ala_plate['time'], y=ala_plate['mCherry'], marker=dict(
        color='blue'), showlegend=False), row=6, col=1)

    fig.add_trace(go.Scatter(x=ala_t_plate['time'], y=ala_t_plate['GFP'], marker=dict(
        color='blue'), showlegend=False), row=5, col=2)
    fig.add_trace(go.Scatter(x=ala_t_plate['time'], y=ala_t_plate['mCherry'], marker=dict(
        color='blue'), showlegend=False), row=6, col=2)


    fig.add_trace(go.Scatter(x=at_ala['sample_time'], y=at_ala['average'],
                            marker=dict(color='#1B9E77'), error_y=dict(type='data',
                            array=at_ala['stdev'].to_list(), visible=True), showlegend=False), row=7, col=1)
    fig.add_trace(go.Scatter(x=at_ala_t['sample_time'], y=at_ala_t['average'],
                            marker=dict(color='#1B9E77'), error_y=dict(type='data',
                            array=at_ala_t['stdev'].to_list(), visible=True), showlegend=False), row=7, col=2)

    fig.add_trace(go.Scatter(x=oa_ala['sample_time'], y=at_ala['average'],
                            marker=dict(color='#D95F02'), error_y=dict(type='data', array=at_ala['stdev'].to_list(), visible=True), showlegend=False), row=7, col=1)
    fig.add_trace(go.Scatter(x=oa_ala_t['sample_time'], y=at_ala_t['average'],
                            marker=dict(color='#D95F02'), error_y=dict(type='data', array=at_ala_t['stdev'].to_list(), visible=True), showlegend=False), row=7, col=2)


    fig['layout']['yaxis13']['type'] = 'log'
    fig['layout']['yaxis14']['type'] = 'log'
    fig['layout']['xaxis13']['title'] = 'Time [h]'
    fig['layout']['xaxis14']['title'] = 'Time [h]'

    return fig

fig = plate_reader()
fig.show()