import plotly.express as px
import pandas as pd
from chibio_parser import *
from os.path import join
from alle_effect import *
from analysis import competition as box_comp
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from citreate_params import get_r
from scipy.integrate import odeint


colors = {'ct': '#7570B3',
          'oa': '#D95F02',
         }

abb = {'ct': '<i>C. testosteroni</i>',
       'oa': '<i>B. anthropi</i>',
       '<i>C. testosteroni</i>': 'ct',
       '<i>B. anthropi</i>': 'oa',
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
                  facet_col='reactor', width=w * 1.2, height=h/1.5)
    titles = [abb['oa'], abb['ct']]
    for i, d in enumerate(fig['data']):
        d['name'] = titles[i]
        d['line']['color'] = colors[abb[titles[i]]]
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    fig.add_vline(
        x=141, col=1, annotation_text='Added Thiamine <br>to feed c=1.48 \u03BCM', opacity=0.7)
    fig.add_vline(x=0, annotation_text='Batch <br>culture', opacity=.7)
    fig.add_vline(x=67, annotation_text='Start <br>dilution', opacity=.7)
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
    fig.add_hline(y=0.1185, line_dash='dash', line_width=2, line_color='black')
    fig = style_plot(fig)
    fig.write_image(join('mid_thesis_plots', '2c.svg'))
    return fig


def fig_2d():
    fig = simulation(0.1, 500)[0]
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


def fig_3():
    fig = make_subplots(rows=2, cols=3, vertical_spacing=0.2, horizontal_spacing=0.12, row_heights=[
                        0.3, 0.7], subplot_titles=['D = 0.04', 'D = 0.15', 'D = 0.085'])
    data_dir = join('..', 'data', 'citrate_thiamine_merged',
                    'dfs', 'M4', '0.04')
    od = pd.read_csv(
        join(data_dir, 'od.csv'))
    fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                             marker=dict(color='#838285'), name='OD', mode="lines", showlegend=False), row=1, col=1)
    for s, color in colors.items():
        cfus = pd.read_csv(join(data_dir, s + '.csv'))
        fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                                 y=cfus['average'], marker=dict(color=color), name=s, showlegend=False), row=2, col=1)

    fig['layout']['width'], fig['layout']['height'] = 2 * w, 2 * h

    data_dir = join('..', 'data', 'citrate_thiamine_merged',
                    'dfs', 'M4', '0.15')
    od = pd.read_csv(
        join(data_dir, 'od.csv'))
    fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                             marker=dict(color='#838285'), name='OD', mode="lines", showlegend=False), row=1, col=2)
    for s, color in colors.items():
        cfus = pd.read_csv(join(data_dir, s + '.csv'))
        fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                                 y=cfus['average'], marker=dict(color=color), name=s, showlegend=False), row=2, col=2)

    data_dir = join('..', 'data', 'citrate_thiamine_merged',
                    'dfs', 'M4', '0.085')
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

    xaxes = ['xaxis' + j for j in ['', '2', '3', '4', '5', '6']]
    yaxes = ['yaxis' + j for j in ['', '2', '3', '4', '5', '6']]

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
    # fig.show()


def supp_fig1():
    fig = make_subplots(rows=2, cols=3, vertical_spacing=0.2, horizontal_spacing=0.12, row_heights=[
        0.3, 0.7], subplot_titles=['D = 0.04', 'D = 0.15', 'D = 0.085'])
    data_dir = join('..', 'data', 'citrate_thiamine_merged',
                    'dfs', 'M5', '0.04')
    od = pd.read_csv(
        join(data_dir, 'od.csv'))
    fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                             marker=dict(color='#838285'), name='OD', mode="lines", showlegend=False), row=1, col=1)
    for s, color in colors.items():
        cfus = pd.read_csv(join(data_dir, s + '.csv'))
        fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                                 y=cfus['average'], marker=dict(color=color), name=s, showlegend=False), row=2, col=1)

    fig['layout']['width'], fig['layout']['height'] = 2 * w, 2 * h

    data_dir = join('..', 'data', 'citrate_thiamine_merged',
                    'dfs', 'M5', '0.15')
    od = pd.read_csv(
        join(data_dir, 'od.csv'))
    fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                             marker=dict(color='#838285'), name='OD', mode="lines", showlegend=False), row=1, col=2)
    for s, color in colors.items():
        cfus = pd.read_csv(join(data_dir, s + '.csv'))
        fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                                 y=cfus['average'], marker=dict(color=color), name=s, showlegend=False), row=2, col=2)

    data_dir = join('..', 'data', 'citrate_thiamine_merged',
                    'dfs', 'M5', '0.085')
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

    xaxes = ['xaxis' + j for j in ['', '2', '3', '4', '5', '6']]
    yaxes = ['yaxis' + j for j in ['', '2', '3', '4', '5', '6']]

    for y in yaxes[:3]:
        fig['layout'][y]['range'] = [0, 0.7]

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
    fig.write_image(join('mid_thesis_plots', 'supp_fig_1.svg'))


def supp_fig2a():
    figs, df = box_comp()
    df = df[(df['species'] == 'ct') & (
        (df['cond'] == 'mono_thiamine') | (df['cond'] == 'mono'))]
    df = df[(df['D'] == 0.04) | (df['D'] == 0.15)]
    fig = px.box(df, x='D', y='CFUs', log_y=True, title='ct_ds', color='species', height=250, width=600, category_orders={
        'D': [0.04,  0.15]}, points='all', facet_col='cond')
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.update_layout(boxgroupgap=0, boxgap=0.2, boxmode='overlay')
    fig.update_xaxes(type='category')
    titles = ['1.48 \u03BCM supplied thiamine', 'No thiamine', '', '']
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    for d in fig['data']:
        d['name'] = 'ct'
        d['line']['color'] = colors[d['name']]
        d['marker']['color'] = colors[d['name']]
    fig['layout']['width'], fig['layout']['height'] = w, h
    fig.update_xaxes(title='Dilution rate 1/h')
    fig.update_yaxes(title='CFUs/mL')
    fig.update_layout(showlegend=False,boxgroupgap=0.5, boxgap=0.3)
    fig = style_plot(fig)
    for d in fig['data']:
        d['marker']['size'] = 4
        d['line']['width'] = 2
    fig.write_image(join('mid_thesis_plots', 'supp_fig_2.svg'))

def supp_fig_2_b():
    df = get_r()[2]

    def model(y, t, q, r, K, D):
        R = y[0]
        N1 = y[1]
        dR = D * M - D * R - N1 / q * r * R / \
            (R + K)
        dN1 = r * R / (K + R) * N1 - D * N1
        return [dR, dN1]

    def model_Oa(y, t, q, r, K, D):
        # Only coded for batch cutlure
        R = y[0]
        T = y[1]
        N1 = y[2]
        J = r * R / (R + K) * T / (T + 0.00074)
        dR = D * M - D * R - J * N1 / q
        dT = -D * T - J * N1 / 831.756756757
        dN1 = J * N1 - D * N1
        return [dR, dT, dN1]

    xs = df[df['species'] == 'ct']['time'].to_list()
    N0 = df[df['species'] == 'ct']['OD'].to_list()[0]
    M = 10
    y = odeint(model, [M, N0], xs, args=(0.065, 0.24, 7, 0))
    ct = pd.DataFrame(columns=['OD', 'time', 'species'])
    ct['OD'], ct['time'], ct['species'] = y[:, 1], xs, 'ct_modelled'

    xs = df[df['species'] == 'oa']['time'].to_list()
    N0 = df[df['species'] == 'oa']['OD'].to_list()[0]
    M = 10
    y = odeint(model_Oa, [M, 0.00148, N0], xs, args=(0.11, 0.43, 7, 0))
    oa = pd.DataFrame(columns=['OD', 'time', 'species'])
    oa['OD'], oa['time'], oa['species'] = y[:, 2], xs, 'oa_modelled'
    out = pd.concat([df, ct, oa])
    out = out.sort_values(by=['species', 'time'])
    fig = px.line(out, x='time', y='OD', color='species')
    l_names = {'ct': abb['ct'],
                'oa': abb['oa'],
                'ct_modelled': abb['ct'] + ' modeled',
                'oa_modelled': abb['oa'] + ' modeled'}
    for d in fig['data']:
        if (d['name'] == 'ct_modelled') | (d['name'] == 'oa_modelled'):
            d['line']['dash'] = 'dash'
        d['line']['color'] = colors[d['name'][:2]]
        d['name'] = l_names[d['name']]
    fig['layout']['width'], fig['layout']['height'] = w, h
    fig.update_xaxes(title='Time [h]')
    fig.update_yaxes(title='OD')
    fig['layout']['legend']['title']['text'] = 'Species'
    fig = style_plot(fig, marker_size=2)
    fig.write_image(join('mid_thesis_plots', 'supp_fig_2b.svg'))
    fig.show()


def plot_all():
    fig_1()
    fig_2a()
    fig_2b()
    fig_2c()
    fig_2d()
    fig_3()
    supp_fig1()
    supp_fig2a()
    supp_fig_2_b()
