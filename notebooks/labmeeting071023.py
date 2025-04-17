from chibio_parser import *
from plotting import plot_od
from os.path import join
import plotly.express as px
from plotting import *
from os import symlink
from sympy import preview
import numpy as np
from cross_feeding import thiamine_cross_feeding
from analysis import plot as exp_plot
from analysis import competition as box_comp
from analysis import competition_line


def update_labels(fig, x_label, y_label, title, size=[250, 400], width=6, style_colors=True, font_size=14, reactor_tiles=False):
    fig.for_each_xaxis(lambda axis: axis.title.update(text=None))
    fig.for_each_yaxis(lambda axis: axis.title.update(text=None))
    fig.for_each_xaxis(lambda axis: axis.title.update(
        font=dict(size=font_size, color='black')))
    fig.for_each_yaxis(lambda axis: axis.title.update(
        font=dict(size=font_size, color='black')))
    lw = 1
    color = '#999999'
    line_color = '#372f55'
    fig.update_layout(yaxis_title=y_label, xaxis_title=x_label, title=title,
                      height=size[0], width=size[1], paper_bgcolor=color, plot_bgcolor=color, font_color='black')
    fig.update_xaxes(showline=True, linewidth=lw,
                     linecolor='black', gridwidth=lw, gridcolor='black', zerolinecolor='black', zerolinewidth=lw)
    fig.update_yaxes(showline=True, linewidth=lw,
                     linecolor='black', gridwidth=lw, gridcolor='black', zerolinecolor='black', zerolinewidth=lw)
    """for anno in fig['layout']['annotations']:
        anno['text'] = ''"""
    # fig.update_layout(title=None)
    fig.update_layout(
        margin=dict(l=10, r=10, t=50, b=10),
    )
    fig.update_layout(font={'size': font_size, 'color': 'black'})
    fig.update_traces(line={'width': width})
    if style_colors:
        fig.update_traces(line={'color': line_color})
    fig['layout']['legend']['title']['text'] = ''
    if reactor_tiles:
        titles = ['Reactor '+str(n) for n in range(1, 5)]
        for counter, annotation in enumerate(fig['layout']['annotations']):
            annotation['text'] = titles[counter]
    for a in fig['layout']['annotations']:
        a['font']['size'] = font_size
        a['font']['color'] = 'black'
    fig['layout']['title']['font']['size'] = font_size
    fig['layout']['title']['font']['color'] = 'black'
    fig['layout']['legend']['title']['font']['size'] = font_size
    fig['layout']['legend']['title']['font']['color'] = 'black'

    return fig


def species_colors(fig):
    colors = {'Ct': '#8872cd',
              'Oa': '#e27e50',
              'ct': '#8872cd',
              'oa': '#e27e50',
              'ct, co': '#8872cd',
              'oa, co': '#e27e50',
              'ct, mono': '#8872cd',
              'oa, mono': '#e27e50'}
    for d in fig['data']:
        d['line']['color'] = colors[d['name']]
        d['marker']['color'] = colors[d['name']]
    return fig


def write_e(expr, f):
    preview(expr, viewer='file', filename=join('figures_labmeeting071223',
            'equations', f), dvioptions=['-D', '600', '-bg', 'Transparent'])


def intro_ct_oa():
    df, o = chibio_parser('ct_oa_citrate_mono')
    df = df[df['sensor'] == 'od_measured']
    df, fig = plot_od('ct_oa_citrate_mono', df=df)
    fig['data'][0]['name'] = 'Oa'
    fig['data'][1]['name'] = 'Ct'
    fig = species_colors(fig)
    fig = update_labels(fig, 'time [h]', 'OD',
                        None, style_colors=False, size=[250, 500])
    fig.update_layout(showlegend=False)
    fig.add_vline(
        x=141, col=1, opacity=1)
    fig.write_image(join('figures_labmeeting071223', 'introduction_ct_oa.svg'))


def batch():
    ct = pd.read_csv(join('..', 'data', 'mono_ss_report', 'ss_17_3_ct.csv'))[
        ['B3']]
    ct['Time'] = np.arange(0, 72.1, 1/6)
    ct['species'] = 'Ct'
    ct.columns = ['OD', 'time', 'species']
    oa = pd.read_csv(join('..', 'data', 'mono_ss_report', 'oa_03_19_curves.csv'))[
        ['Time', 'Citric acid']]
    oa['species'] = 'Oa'
    oa.columns = ['time', 'OD', 'species']
    oa['time'] = np.arange(0, 44.6, 1/6)
    out = pd.concat([ct, oa])
    out = out[out['time'] <= 40]
    fig = px.line(out, x='time', y='OD', color='species')
    fig = species_colors(fig)
    fig = update_labels(fig, 'time [h]', 'OD', None,
                        style_colors=False, width=4, size=[250, 350])
    fig.write_image(join('figures_labmeeting071223', 'batch_ct_oa.svg'))
    return fig


def equations():
    K = '$K = \\frac{R(\\mu_{max} - D)}{D}$'
    write_e(K, 'Km.svg')
    S = '$R = \\frac{-(N - YR_{feed})}{Y}$'
    write_e(S, 'S.svg')
    J_c_ct = '$J_{C,Ct}^{upt} = \\mu_{Ct}\\frac{C}{K_{C,Ct}+C}$'
    write_e(J_c_ct, 'j_c_ct.svg')
    J_c_oa = '$J_{C,Oa}^{upt} = \\mu_{Oa}\\frac{C}{K_{C,Oa}+C}$'
    write_e(J_c_oa, 'j_c_oa.svg')
    J_t_oa = '$J_{T,Oa}^{upt} = \\mu_{Oa}\\frac{T}{K_{T,Oa}+T}$'
    write_e(J_t_oa, 'j_t_oa.svg')
    dC = '$\\frac{dC}{dt} = D(M_C-C)-J_{C,Ct}^{upt}\\frac{Ct}{Y_{C,Ct}}-J_{C,Oa}^{upt}\\frac{Oa}{Y_{C,Oa}}$'
    write_e(dC, 'dC.svg')
    dOa = '$\\frac{dT}{dt} = \\alpha J_{C,Ct}^{upt}\\frac{Ct}{Y_{C,Ct}} - J_{T,Oa}^{upt}\\frac{Oa}{Y_{T,Oa}}-DT$'
    write_e(dOa, 'dOa.svg')
    dCt = '$\\frac{dCt}{dt} = Ct(J_{C,Ct}^{upt} - D)$'
    write_e(dCt, 'dCt.svg')
    dT = '$\\frac{dT}{dt} = \\alpha J_{C,Ct}^{upt}\\frac{Ct}{Y_{C,Ct}} - J_{T,Oa}^{upt}\\frac{Oa}{Y_{T,Oa}}-DT$'
    write_e(dT, 'dT.svg')
    dOa = '$\\frac{dOa}{dt} = Oa(min(J_{C,Oa}^{upt},J_{T,Oa}^{upt}) - D)$'
    write_e(dOa, 'dOa.svg')


def cross_feeding():
    f, c, t, y = thiamine_cross_feeding(0.04, range(800))
    f = update_labels(f, 'time [h]', 'CFUs/mL',
                      'D=0.04', style_colors=False, width=3)
    f.write_image(join('figures_labmeeting071223', 'cross_feeding_0.04.svg'))
    f, c, t, y = thiamine_cross_feeding(0.15, range(800))
    f = update_labels(f, 'time [h]', 'CFUs/mL',
                      'D=0.15', style_colors=False, width=3)
    f.write_image(join('figures_labmeeting071223', 'cross_feeding_0.15.svg'))
    f, c, t, y = thiamine_cross_feeding(0.085, range(800))
    f = update_labels(f, 'time [h]', 'CFUs/mL',
                      'D=0.085', style_colors=False, width=3)
    f.write_image(join('figures_labmeeting071223', 'cross_feeding_0.085.svg'))


def chem_exps():
    figs = exp_plot()
    fig = update_labels(figs[0], 'time [h]', 'Population size OD or CFUs/mL',
                        None, size=[450, 1000], style_colors=False, width=3)
    fig.write_image(join('figures_labmeeting071223', 'd_04.svg'))
    fig = update_labels(figs[1], 'time [h]', 'Population size OD or CFUs/mL',
                        None, size=[450, 1000], style_colors=False, width=3)
    fig.write_image(join('figures_labmeeting071223', 'd_15.svg'))
    fig = update_labels(figs[2], 'time [h]', 'Population size OD or CFUs/mL',
                        None, size=[450, 1000], style_colors=False, width=3)
    fig.write_image(join('figures_labmeeting071223', 'd_085.svg'))
    fig = update_labels(figs[3], 'time [h]', 'Population size OD or CFUs/mL',
                        None, size=[450, 1000], style_colors=False, width=3)
    fig.write_image(join('figures_labmeeting071223', 'd_all.svg'))



def interactions():
    figs = competition_line()
    fig = update_labels(figs[0], 'time [h]', 'CFUs/mL', None,
                        size=[300, 600], style_colors=False, width=5)
    fig = species_colors(fig)
    fig.write_image(join('figures_labmeeting071223', 'box_0.04.svg'))
    fig = update_labels(figs[1], 'time [h]', 'CFUs/mL', None,
                        size=[300, 600], style_colors=False, width=5)
    fig = species_colors(fig)
    fig.write_image(join('figures_labmeeting071223', 'box_0.15.svg'))
    fig = update_labels(figs[2], 'time [h]', 'CFUs/mL', None,
                        size=[300, 600], style_colors=False, width=5)
    fig = update_labels(figs[3], 'time [h]', 'CFUs/mL', None,
                        size=[300, 600], style_colors=False, width=5)
    fig.write_image(join('figures_labmeeting071223', 'box_all.svg'))
    fig = species_colors(fig)
    fig.write_image(join('figures_labmeeting071223', 'box_0.085.svg'))
    figs = box_comp()
    fig = update_labels(figs[3],None,'CFUs/mL',None,size=[250,600],style_colors=False,width=3)
    fig = species_colors(fig)
    fig.write_image(join('figures_labmeeting071223','box_ds.svg'))
    fig = update_labels(figs[4],None,'CFUs/mL',None,size=[250,400],style_colors=False,width=3)
    fig.write_image(join('figures_labmeeting071223','box_total.svg'))

interactions()