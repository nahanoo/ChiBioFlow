from plotting import *
import plotly.express as px
from chibio_parser import cfu_parser
from os.path import join
from community_caller import update_labels
from citrate import constant_thiamine


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
    fig.write_image(join('report_figures', 'od_comm.png'))
    return fig


def cfus_comm():
    df = cfu_parser('community_ss_report')[0]
    titles = {'M1': 'Glucose', 'M2': 'Citric acid',
              'M3': 'Glucose + Citric acid'}
    for i, r in enumerate(df['reactor']):
        df.at[i, 'reactor'] = titles[r]
    cfus = px.scatter(df, x='reactor', y='average', color='species',
                      category_orders={'reactor': list(titles.values())}, error_y='stdev', log_y=True)
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
    batch['time'], batch['OD'], batch['condition'] = f['Time'], f['OD'], 'MM after overnight'
    cont = pd.DataFrame(columns=['time', 'OD', 'condition'])
    f = pd.read_csv(join('data', 'mono_ss_report', 'ss_17_3_ct.csv'))[['B3']]
    f['Time'] = np.arange(0, 72.1, 1/6)
    cont['time'], cont['OD'], cont['condition'] = f['Time'], f['B3'], 'MM after continuous'
    out = pd.concat([batch, cont])
    fig = px.line(out, x='time', y='OD', color='condition')
    fig.update_layout(autosize=False,
                      width=400,
                      height=300,)
    fig.update_xaxes(title='Time [h]')
    colors = ['#b3c2f2', '#735cdd']
    for d in fig['data']:
        if d['name'] == 'MM after overnight':
            d['line']['color'] = colors[0]
        else:
            d['line']['color'] = colors[1]
    fig.update_layout(title='<i>C. testosteroni</i> Citric acid')
    fig.write_image(join('report_figures', 'ct_plate_reader.png'))
    return fig


def citrate_0_1():
    df = cfu_parser('citrate_ss')[0]
    for i, r in enumerate(df['reactor']):
        if r == "M2":
            df.at[i, 'reactor'] = 'Community'
        if r == "M0":
            df.at[i, 'reactor'] = '<i>C. testosteroni</i>'
    cfus = px.scatter(df, x='reactor', y='average',
                      color='species', error_y='stdev',
                      category_orders={'reactor': ['<i>C. testosteroni</i>', 'Community']}, log_y=True, range_y=[100E6, 5E9])
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
    titles = ['Community','<i>C. testosteroni</i>']
    fig = plot_od('citrate_ss')[1]
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    for d in fig['data']:
        d['line']['color'] = '#372f55'
    fig.update_layout(showlegend=False, autosize=True,
                      width=600,)
    fig.write_image(join('report_figures', 'citrate_0.1_OD.png'))
    return cfus, rel, fig


def mono(subset=False):
    colors = {'M0': '#2c8c5a',
              'M1': '#e27e50',
              'M2': '#e27e50',
              'M3': '#e27e50',
              'M4': '#8872cd'}
    df, fig = plot_od('mono_ss_report')

    titles = ['Glucose, D=0.15', 'Glucose, D=0.05',
              'Citric acid, D=0.05',
              'Glucose + Citric acid, D=0.05',
              'Citric acid, D=0.05']
    Ms = ['M0', 'M1', 'M2', 'M3', 'M4']
    M_dict = {M: titles[i] for i, M in enumerate(Ms)}
    """for i, r in enumerate(df['reactor']):
        df.at[i, 'reactor'] = Ms[r]"""
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    for d in fig['data']:
        d['line']['color'] = colors[d['name']]
    fig.update_layout(showlegend=False, autosize=True,
                      width=1000,)
    fig.write_image(join('report_figures', 'od_mono.png'))

    df = cfu_parser('mono_ss_report')[0]
    titles = ['Glucose, D=0.15', 'Glucose, D=0.05', 'Citric acid, D=0.05',
              'Glucose + Citric acid, D=0.05', 'Citric acid, D=0.05']
    cfus = px.scatter(df, x='reactor', y='average',
                      color='species', error_y='stdev',
                      log_y=True, category_orders={'reactor': Ms})
    cfus = style_plot(None, cfus, 'cfus')
    cfus.update_layout(width=550, xaxis=dict(tickmode='array',
                                             tickvals=[0, 1, 2, 3, 4],
                                             ticktext=titles))
    cfus.update_xaxes(title=None)
    cfus.write_image(join('report_figures', 'cfus_mono.png'))
    return fig


def oa_batch_cont():
    f = pd.read_csv(join('data', 'plate_reader_resources', 'data.csv'))
    B6 = f[f['well'] == 'B6']
    gluc = pd.DataFrame(columns=['time', 'OD', 'condition', 'trans'])
    gluc['time'], gluc['OD'], gluc['condition'], gluc['trans'] = B6['Time'], B6['OD'], 'Glucose', 'MM after overnight'
    cit = pd.DataFrame(columns=['time', 'OD', 'condition', 'trans'])
    C6 = f[f['well'] == 'C6']
    cit['time'], cit['OD'], cit['condition'], cit['trans'] = C6['Time'], C6['OD'], 'Citric acid', 'MM after overnight'
    D6 = f[f['well'] == 'D6']
    gluc_cit = pd.DataFrame(columns=['time', 'OD', 'condition', 'trans'])
    gluc_cit['time'], gluc_cit['OD'], gluc_cit['condition'], gluc_cit['trans'] = D6['Time'], D6['OD'], 'Glucose + Citric acid', 'MM after overnight'
    cont_cit = pd.DataFrame(columns=['time', 'OD', 'condition', 'trans'])
    f = pd.read_csv(join('data', 'mono_ss_report', 'oa_03_19_curves.csv'))
    f['Time'] = np.arange(0, 44.6, 1/6)
    c = "MM + 1.8" + u"\u03BC"+"M Thiamine after continuous"
    cont_cit['time'], cont_cit['OD'], cont_cit['condition'], cont_cit['trans'] = f['Time'], f['Citric acid'], 'Citric acid', c
    cont_gluc = pd.DataFrame(columns=['time', 'OD', 'condition', 'trans'])
    cont_gluc['time'], cont_gluc['OD'], cont_gluc['condition'], cont_gluc['trans'] = f['Time'], f['Glucose'], 'Glucose', c
    cont_cit_gluc = pd.DataFrame(columns=['time', 'OD', 'condition', 'trans'])
    cont_cit_gluc['time'], cont_cit_gluc['OD'], cont_cit_gluc['condition'], cont_cit_gluc[
        'trans'] = f['Time'], f['Citric acid + Glucose'], 'Glucose + Citric acid', c

    out = pd.concat([gluc, cit, gluc_cit, cont_gluc, cont_cit, cont_cit_gluc])
    fig = px.line(out, x='time', y='OD', facet_col='condition', color='trans')
    titles = ['Glucose', 'Citric acid', 'Glucose + Citric acid']
    fig.update_layout(autosize=False,
                      width=700,
                      height=300,)
    fig.update_xaxes(title='Time [h]')
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    colors = ['#b3c2f2', '#735cdd']
    for d in fig['data']:
        if d['name'] == 'MM after overnight':
            d['line']['color'] = colors[0]
        else:
            d['line']['color'] = colors[1]
    fig.update_layout(title='<i>O. anthropi</i>')
    fig.write_image(join('report_figures', 'oa_plate_reader.png'))
    return fig


def ct_oa():
    colors = {'<i>C. testosteroni</i>': '#8872cd',
              '<i>O. anthropi</i>': '#e27e50'}
    df, fig = plot_od('ct_oa_citrate_mono')

    titles = ['Citric acid, D=0.05',
              'Citric acid, D=0.05']
    Ms = ['M0', 'M1', 'M2', 'M3', 'M4']
    names = ['<i>O. anthropi</i>','<i>C. testosteroni</i>']
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    for i, d in enumerate(fig['data']):
        d['name'] = names[i]
        d['line']['color'] = colors[d['name']]
    fig.add_vline(
        x=141, col=1, annotation_text='Added Thiamine <br>to feed c=1.88 \u03BCM', opacity=1)
    fig.add_vline(x=0, annotation_text='Batch <br>culture', opacity=1)
    fig.add_vline(x=67, annotation_text='Start <br>dilution', opacity=1)
    fig['layout']['legend']['title']['text'] = 'species'
    return fig


def od_citric_acid_comm():
    titles = ['D=0.1', 'D=0.15']
    fig = plot_od('citric_acid_comm')[1]
    for i, d in enumerate(fig['data']):
        d['line']['color'] = '#372f55'
    for i, a in enumerate(fig['layout']['annotations']):
        a['text'] = titles[i]
    fig.update_layout(showlegend=False)
    fig.add_vline(x=45, annotation_text='Start <br>dilution', opacity=0.5)
    fig.show()
    return fig

def at_monot():
    fig = plot_od('at_mono')[1]
    fig['layout']['annotations'][0]['text'] = ''
    fig = update_labels(fig,'Time [h]','OD','<i>A. tumefaciens</i> in Glucose',size=[250,400],style_colors=True,width=3)
    fig.update_layout(showlegend=False)
    fig.show()

dfs = []
df015 = cfu_parser('community_ss_report')[0]
df015.insert(len(df015.columns), 'condition', 'Community')
df015.insert(len(df015.columns), 'D', 0.15)
df015 = df015.loc[[0, 6]]
df015['reactor'] = 'D=0.15'
dfs.append(df015)
df01 = cfu_parser('citrate_ss')[0]
df01.insert(len(df01.columns), 'condition', 'Community')
df01.insert(len(df01.columns), 'D', 0.1)
df01 = df01.loc[[0, 2]]
df01['reactor'] = 'D=0.1'
dfs.append(df01)
out = pd.concat(dfs)
fig = px.scatter(out,x='reactor',y='composition',color='species',category_orders={'reactor':['D=0.1','D=0.15']})
fig = style_plot(None,fig,'cfus')
