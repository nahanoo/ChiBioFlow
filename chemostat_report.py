from plotting import *
import plotly.express as px
from chibio_parser import cfu_parser
from os.path import join


def od_comm():
    titles = {'M1': 'Glucose', 'M2': 'Citric acid',
              'M3': 'Glucose + Citric acid'}
    fig = plot_od('community_ss_report')[1]
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = list(titles.values())[i]
    for d in fig['data']:
        d['line']['color'] = '#372f55'
    fig.update_layout(showlegend=False, autosize=False,
                      width=600,
                      height=300,)
    return fig


def cfus_comm():
    df = cfu_parser('community_ss_report')[0]
    titles = {'M1': 'Glucose', 'M2': 'Citric acid',
              'M3': 'Glucose + Citric acid'}
    for i, r in enumerate(df['reactor']):
        df.at[i, 'reactor'] = titles[r]
    cfus = px.scatter(df, x='reactor', y='average', color='species',
                      category_orders={'reactor': list(titles.values())}, error_y='stdev')
    cfus = style_plot(None, cfus, 'cfus')
    cfus.update_layout(width=350)
    cfus.write_image(join('report_figures', 'community_abs_abun.png'))

    rel = px.scatter(df, x='reactor', y='composition', color='species',
                     category_orders={'reactor': list(titles.values())})
    rel = style_plot(None, rel, 'cfus')
    rel.update_layout(width=350)
    rel.update_yaxes(title='Relative abundance')
    rel.write_image(join('report_figures', 'community_rel_abun.png'))
    return cfus, rel


def od_mono():
    titles = ['At Glucose D=0.15', 'Oa Glucose D=0.05',
              'Oa Citric acid D=0.05', 'Oa Glucose + Citric acid D=0.05', 'Ct Citric acid D=0.05']
    fig = plot_od('mono_ss_report')[1]
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    for d in fig['data']:
        d['line']['color'] = '#372f55'
    fig.update_layout(showlegend=False, autosize=True,
                      height=300,)
    return fig


def ct_batch_cont():
    batch = pd.DataFrame(columns=['time', 'OD', 'condition'])
    f = pd.read_csv(join('data', 'plate_reader_resources', 'data.csv'))
    f = f[f['well'] == 'F4']
    batch['time'], batch['OD'], batch['condition'] = f['Time'], f['OD'], 'Batch'
    cont = pd.DataFrame(columns=['time', 'OD', 'condition'])
    f = pd.read_csv(join('data', 'mono_ss_report', 'ss_17_3_ct.csv'))[['B3']]
    f['Time'] = np.arange(0, 72.1, 1/6)
    cont['time'], cont['OD'], cont['condition'] = f['Time'], f['B3'], 'Continuous'
    out = pd.concat([batch, cont])
    fig = px.line(out, x='time', y='OD', color='condition')
    fig.update_layout(autosize=False,
                      width=400,
                      height=300,)
    fig.update_xaxes(title='Time [h]')
    return fig


def citrate_0_1():
    df = cfu_parser('citrate_ss')[0]
    for i, r in enumerate(df['reactor']):
        if r == "M2":
            df.at[i, 'reactor'] = 'Community'
        if r == "M0":
            df.at[i, 'reactor'] = '<i>C. testosteroni</i>'
    cfus = px.scatter(df, x='reactor', y='average',
                      color='species', error_y='stdev', category_orders={'reactor': ['<i>C. testosteroni</i>', 'Community']})
    cfus = style_plot(None, cfus, 'cfus')
    cfus.update_layout(width=350, title='Citric acid')
    cfus.update_xaxes(title=None)
    cfus.write_image(join('report_figures', 'citrate_0.1_abs_abun.png'))
    df = df[df['reactor'] == 'Community']
    rel = px.scatter(df, x='reactor', y='composition', color='species')
    rel = style_plot(None, rel, 'cfus')
    rel.update_yaxes(title='Relative abundance')
    rel.update_xaxes(title=None)
    for d in rel['data']:
        d['x'][0] = 'Community'
    rel.update_layout(width=350, title='Citric acid')

    rel.write_image(join('report_figures', 'citrate_0.1_rel_abun.png'))
    titles = ['<i>C. testosteroni</i>','Community']
    fig = plot_od('citrate_ss')[1]
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    for d in fig['data']:
        d['line']['color'] = '#372f55'
    fig.update_layout(showlegend=False, autosize=True,
                      width=600,)
    fig.write_image(join('report_figures', 'citrate_0.1_OD.png'))
    return cfus, rel, fig