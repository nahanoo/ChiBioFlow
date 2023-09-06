import plotly.express as px
import pandas as pd
from chibio_parser import *
from os.path import join
from alle_effect import *
from plotly.subplots import make_subplots
import plotly.graph_objects as go

colors = {'ct': '#7570B3',
          'oa': '#D95F02'}

abb = {'ct': '<i>C. testosteroni</i>',
       'oa': '<i>O. anthropi</i>',
       '<i>C. testosteroni</i>': 'ct',
       '<i>O. anthropi</i>': 'oa',
       'OD': 'OD'
       }

w = 400
h = 250


def style_plot(fig, marker_size=3):
    """Style function for figures setting fot size and true black color."""
    for d in fig['data']:
        d['marker']['size'] = marker_size
        d['line']['width'] = marker_size
    # Font size
    j = 10
    fig.update_layout(font={'size': j, 'color': 'black'})
    for a in fig['layout']['annotations']:
        a['font']['size'] = j
        a['font']['color'] = 'black'
    fig['layout']['title']['font']['size'] = j
    fig['layout']['title']['font']['color'] = 'black'
    fig['layout']['legend']['title']['font']['size'] = j
    fig['layout']['legend']['title']['font']['color'] = 'black'
    fig.for_each_xaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    fig.for_each_yaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    fig.update_layout(
        margin=dict(l=60, r=20, t=20, b=20),
    )
    fig.update_yaxes(title_standoff=10)
    fig.update_xaxes(title_standoff=10)
    return fig


def fig_1():
    df, o = chibio_parser('ct_oa_citrate_mono')
    df = df[df['sensor'] == 'od_measured']
    fig = px.line(df, x='exp_time', y='measurement',
                  facet_col='reactor', width=w*2, height=h)
    titles = [abb['oa'], abb['ct']]
    for i, d in enumerate(fig['data']):
        d['name'] = titles[i]
        d['line']['color'] = colors[abb[titles[i]]]
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    fig.add_vline(
        x=141, col=1, annotation_text='Added Thiamine <br>to feed c=1.48 \u03BCM', opacity=1)
    fig.add_vline(x=0, annotation_text='Batch <br>culture', opacity=1)
    fig.add_vline(x=67, annotation_text='Start <br>dilution', opacity=1)
    fig['layout']['legend']['title']['text'] = 'species'
    fig.update_xaxes(title='Time [h]')
    fig.update_yaxes(title='OD')
    fig = style_plot(fig, marker_size=1.5)
    fig.update_xaxes(matches=None)
    fig.write_image(join('mid_thesis_plots', '1_a.svg'))
    return fig


def fig_2a():
    fig = alle_effect()
    fig['layout']['width'], fig['layout']['height'] = w, h
    for d in fig['data']:
        d['line']['color'] = colors[d['name']]
        d['name'] = abb[d['name']]
    fig.update_xaxes(title='Community OD', dtick=0.2, tickmode='linear')
    fig.update_yaxes(
        title='Per capita growth rate [1/h]', dtick=0.02, tickmode='linear')
    fig['layout']['legend']['title']['text'] = 'Species'
    fig = style_plot(fig)
    fig.write_image(join('mid_thesis_plots', '2a.svg'))


def fig_2b():
    fig = composition()
    fig['layout']['width'], fig['layout']['height'] = w, h
    for d in fig['data']:
        d['line']['color'] = colors[d['name']]
        d['name'] = abb[d['name']]
    fig.update_xaxes(
        title='Dilution rate [1/h]', dtick=0.02, tickmode='linear')
    fig.update_yaxes(title='Community composition',
                     dtick=0.2, tickmode='linear')
    fig['layout']['legend']['title']['text'] = 'Species'
    fig = style_plot(fig)
    fig.write_image(join('mid_thesis_plots', '2b.svg'))


def fig_2c():
    fig = alle_threshold()
    fig['layout']['width'], fig['layout']['height'] = w, h
    for d in fig['data']:
        d['line']['color'] = colors[d['name']]
        d['name'] = abb[d['name']]
    fig.update_xaxes(
        title='Dilution rate [1/h]', dtick=0.02, tickmode='linear')
    fig.update_yaxes(
        title='Per capita growth rate [1/h]', dtick=0.02, tickmode='linear')
    fig['layout']['legend']['title']['text'] = 'Species'
    fig = style_plot(fig)
    fig.write_image(join('mid_thesis_plots', '2c.svg'))


def fig_2d():
    fig = simulation(0.09, 500)
    fig['layout']['width'], fig['layout']['height'] = w, h
    for d in fig['data']:
        d['line']['color'] = colors[d['name']]
        d['name'] = abb[d['name']]
    fig.update_xaxes(title='Time [h]', dtick=100, tickmode='linear')
    fig.update_yaxes(title='OD', dtick=0.1, tickmode='linear')
    fig['layout']['legend']['title']['text'] = 'Species'
    fig = style_plot(fig)
    fig.write_image(join('mid_thesis_plots', '2da.svg'))

    fig = simulation_mac(0.09, 500)
    fig['layout']['width'], fig['layout']['height'] = w, h
    for d in fig['data']:
        d['line']['color'] = colors[d['name']]
        d['name'] = abb[d['name']]
    fig.update_xaxes(title='Time [h]', dtick=100, tickmode='linear')
    fig.update_yaxes(title='OD', dtick=0.5, tickmode='linear')
    fig['layout']['legend']['title']['text'] = 'Species'
    fig = style_plot(fig)
    fig.write_image(join('mid_thesis_plots', '2db.svg'))
    return fig


def plotter_2():
    fig_2a()
    fig_2b()
    fig_2c()
    fig_2d()


fig = make_subplots(rows=2, cols=3,vertical_spacing=0.2,row_heights=[0.3,0.7],subplot_titles=['D = 0.04','D = 0.15','D = 0.085'])
data_dir = join('..', 'data', 'citrate_thiamine_merged', 'dfs', 'M4', '0.04')
od = pd.read_csv(
    join(data_dir, 'od.csv'))
fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                         marker=dict(color='#838285'), name='OD', mode="lines",showlegend=False), row=1, col=1)
for s, color in colors.items():
    cfus = pd.read_csv(join(data_dir, s + '.csv'))
    fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                             y=cfus['average'], marker=dict(color=color), name=s,showlegend=False), row=2, col=1)

fig['layout']['width'], fig['layout']['height'] = 2* w, 2 * h

data_dir = join('..', 'data', 'citrate_thiamine_merged', 'dfs', 'M4', '0.15')
od = pd.read_csv(
    join(data_dir, 'od.csv'))
fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                         marker=dict(color='#838285'), name='OD', mode="lines",showlegend=False), row=1, col=2)
for s, color in colors.items():
    cfus = pd.read_csv(join(data_dir, s + '.csv'))
    fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                             y=cfus['average'], marker=dict(color=color), name=s,showlegend=False), row=2, col=2)


data_dir = join('..', 'data', 'citrate_thiamine_merged', 'dfs', 'M4', '0.085')
od = pd.read_csv(
    join(data_dir, 'od.csv'))
fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                         marker=dict(color='#838285'), name='OD', mode="lines"), row=1, col=3)
for s, color in colors.items():
    cfus = pd.read_csv(join(data_dir, s + '.csv'))
    fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                             y=cfus['average'], marker=dict(color=color), name=s), row=2, col=3)

for d in fig['data']:
    d['name'] = abb[d['name']]

xaxes = ['xaxis' + j for j in ['','2','3','4','5','6']]
yaxes = ['yaxis' + j for j in ['','2','3','4','5','6']]


for y in yaxes[:3]:
    fig['layout'][y]['range'] = [0, 0.6]

for y in yaxes[3:]:
    fig['layout'][y]['type'] = 'log'
    fig['layout'][y]['range'] = [6, 10]

fig['layout']['yaxis']['title']['text'] = 'OD'
fig['layout']['yaxis4']['title']['text'] = 'CFUs/mL'


for x in xaxes:
    fig['layout'][x]['tickmode'] = 'linear'
    fig['layout'][x]['dtick'] = 30

fig['layout']['xaxis']['range'] = [110, 215]
fig['layout']['xaxis4']['range'] = [110, 215]
fig['layout']['xaxis4']['title']['text'] = 'Time [h]'


fig['layout']['xaxis2']['range'] = [238, 381]
fig['layout']['xaxis5']['range'] = [238, 381]
fig['layout']['xaxis5']['title']['text'] = 'Time [h]'


fig['layout']['xaxis3']['range'] = [450, 619]
fig['layout']['xaxis6']['range'] = [450, 619]
fig['layout']['xaxis6']['title']['text'] = 'Time [h]'

fig = style_plot(fig)
fig['layout']['width'], fig['layout']['height'] = w*2, h
fig.write_image(join('mid_thesis_plots', '3.svg'))
fig.show()
"""
fig['layout']['yaxis3']['range'] = [0, 0.6]
fig['layout']['yaxis6']['type'] = 'log'
fig['layout']['yaxis6']['range'] = [6, 10]
fig['layout']['xaxis6']['title']['text'] = 'Time [h]'
fig['layout']['xaxis3']['range'] = [450, 619]
fig['layout']['xaxis6']['range'] = [450, 619]

fig['layout']['yaxis2']['range'] = [0, 0.6]
fig['layout']['yaxis5']['type'] = 'log'
fig['layout']['yaxis5']['range'] = [6, 10]
fig['layout']['xaxis5']['title']['text'] = 'Time [h]'
fig['layout']['xaxis5']['range'] = [238, 381]
fig['layout']['xaxis3']['range'] = [238, 381]

fig['layout']['yaxis']['range'] = [0, 0.6]
fig['layout']['yaxis']['title']['text'] = 'OD'
fig['layout']['yaxis4']['type'] = 'log'
fig['layout']['yaxis4']['range'] = [6, 10]
fig['layout']['yaxis4']['title']['text'] = 'CFUs/mL'
fig['layout']['xaxis4']['title']['text'] = 'Time [h]'
fig['layout']['xaxis']['range'] = [110, 215]
fig['layout']['xaxis4']['range'] = [110, 215]"""