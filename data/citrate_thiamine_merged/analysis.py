import plotly.express as px
import pandas as pd
from glob import glob
from os.path import join, split, exists
from chibio_parser import *
from plotting import *
from os import mkdir
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def add_start_time():
    for f in glob(join('..', 'citrate_thiamine', 'M*', '*.csv')):
        df = pd.read_csv(f)
        out = pd.DataFrame(columns=df.columns)
        for i in range(24 * 60 + 18):
            row = [None] * len(df.columns)
            row[0] = i * 60
            out.loc[i] = row
        if split(f)[0][-2:] in ['M1', 'M3', 'M4', 'M5']:
            df['exp_time'] += (24 * 60 + 18) * 60
            out = pd.concat([out, df])
        else:
            out = df
        out_f = join(split(f)[0], 'mutated_' + split(f)[1])
        out.to_csv(out_f, index=False)


def merge():
    ct = []
    for f in glob(join('..', 'citrate_thiamine', 'M*', 'mutated_*.csv')):
        ct.append(f)
    ct = sorted(ct)

    ct_r = []
    for f in glob(join('..', 'citrate_thiamine_restart', 'M*', '*.csv')):
        ct_r.append(f)
    ct_r = sorted(ct_r)

    exp_times = []
    for f in ct:
        t = pd.read_csv(f)[['exp_time']]['exp_time'].to_list()[-1]
        exp_times.append(t)

    max_exp_time = max(exp_times)

    for i, j in zip(ct, ct_r):
        df_i = pd.read_csv(i)
        df_j = pd.read_csv(j)
        df_j['exp_time'] = df_j['exp_time'] + max_exp_time
        out = pd.concat([df_i, df_j])
        reactor = i.split('/')[2]
        f = join(reactor, split(i)[-1])
        out.to_csv(f, index=False)


def competition():
    figs = []
    dfs = []
    out = pd.DataFrame(
        columns=['D', 'cond', 'CFUs', 'sample_time', 'total', 'species'])
    df, o = cfu_parser('citrate_thiamine_merged')
    # df.astype({'sample_time':'float'})
    df = df[df['sample_time'] < 213]
    conditions = {'M0': 'mono_thiamine', 'M2': 'mono',
                  'M4': 'comm', 'M5': 'comm_thiamine'}
    for r, s, c, t, a in zip(df['reactor'], df['species'], df['average'], df['sample_time'], df['total']):
        if (r in conditions.keys()):
            out.loc[len(out)] = [0.04, conditions[r], c, t, a, s]
    out.loc[len(out)] = [0.04, 'mono_thiamine', 1.43E9, 115, None, 'oa']
    out.loc[len(out)] = [0.04, 'mono_thiamine', 1.4E9, 120, None, 'oa']

    out.loc[len(out)] = [0.04, 'mono', 7E8, 115, None, 'oa']
    out.loc[len(out)] = [0.04, 'mono', 6.3E8, 120, None, 'oa']

    fig = px.box(out, x='cond', hover_data=['sample_time'], y='CFUs', log_y=True, title='ct_04', category_orders={
        'cond': ['mono', 'comm', 'mono_thiamine', 'comm_thiamine']}, points='all', height=250, width=400, color='species')
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.update_layout(boxgroupgap=0, boxgap=0.2, boxmode='overlay')
    dfs.append(out)
    figs.append(fig)

    out = pd.DataFrame(columns=['D', 'cond', 'CFUs', 'total', 'species'])
    df, o = cfu_parser('citrate_thiamine_merged')
    # df.astype({'sample_time':'float'})
    df = df[(df['sample_time'] > 213) & (df['sample_time'] < 381)]
    conditions = {'M0': 'mono_thiamine', 'M2': 'mono',
                  'M4': 'comm', 'M5': 'comm_thiamine'}
    for r, s, c, a in zip(df['reactor'], df['species'], df['average'], df['total']):
        if (r in conditions.keys()) & (s == 'ct'):
            out.loc[len(out)] = [0.15, conditions[r], c, a, s]

    fig = px.box(out, x='cond', y='CFUs', log_y=True, title='ct_0.15', color='species', category_orders={
        'cond': ['mono', 'comm', 'mono_thiamine', 'comm_thiamine']}, points='all', height=250, width=400)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    figs.append(fig)
    dfs.append(out)

    out = pd.DataFrame(columns=['D', 'cond', 'CFUs', 'total', 'species'])
    df, o = cfu_parser('citrate_thiamine_merged')
    # df.astype({'sample_time':'float'})
    df = df[df['sample_time'] > 381]
    conditions = {'M0': 'mono_thiamine', 'M2': 'mono',
                  'M4': 'comm', 'M5': 'comm_thiamine'}
    for r, s, c, a in zip(df['reactor'], df['species'], df['average'], df['total']):
        if (r in conditions.keys()) & (s == 'ct'):
            out.loc[len(out)] = [0.085, conditions[r], c, a, s]

    fig = px.box(out, x='cond', y='CFUs', log_y=True, title='ct_0.085', color='species', category_orders={
        'cond': ['mono', 'comm', 'mono_thiamine', 'comm_thiamine']}, points='all', height=250, width=400)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    dfs.append(out)
    figs.append(fig)
    df = pd.concat(dfs)
    ct = df[df['species'] == 'ct']
    fig = px.box(ct, x='D', y='CFUs', log_y=True, title='ct_ds', color='species', height=250, width=600, category_orders={
                 'D': [0.04, 0.085, 0.15]}, points='all', facet_col='cond')
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.update_layout(boxgroupgap=0, boxgap=0.2, boxmode='overlay')
    fig.update_xaxes(type='category')
    figs.append(fig)
    total = df[(df['cond'] == 'comm') | (
        df['cond'] == 'comm_thiamine')].drop_duplicates()
    fig = fig = px.box(total, x='D', y='total', log_y=True, title='total_comm', height=250, width=600, category_orders={
        'D': [0.04, 0.085, 0.15]}, points='all', facet_col='cond')
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.update_layout(boxgroupgap=0, boxgap=0.2, boxmode='overlay')
    fig.update_xaxes(type='category')
    figs.append(fig)
    return figs

def competition_line():
    figs = []
    out = pd.DataFrame(
        columns=['D', 'cond', 'CFUs', 'sample_time', 'total', 'species', 'media', 'mono'])
    df, o = cfu_parser('citrate_thiamine_merged')
    # df.astype({'sample_time':'float'})
    conditions = {'M0': 'mono_thiamine', 'M2': 'mono',
                'M4': 'comm', 'M5': 'comm_thiamine'}
    for r, s, c, t, a in zip(df['reactor'], df['species'], df['average'], df['sample_time'], df['total']):
        if (r in conditions.keys()):
            if 'thiamine' in conditions[r]:
                m = 'Citric acid + Thiamine'
            else:
                m = 'Citric acid'
            if 'mono' in conditions[r]:
                presence = 'mono'
            else:
                presence = 'co'
            out.loc[len(out)] = [0.04, conditions[r], c, t, a, s, m, presence]
    out.loc[len(out)] = [0.04, 'mono_thiamine', 1.43E9, 115,
                        None, 'oa', 'Citric acid + Thiamine', 'mono']
    out.loc[len(out)] = [0.04, 'mono_thiamine', 1.4E9, 120,
                        None, 'oa', 'Citric acid + Thiamine', 'mono']

    out.loc[len(out)] = [0.04, 'mono', 7E8, 115, None, 'oa', 'Citric acid', 'mono']
    out.loc[len(out)] = [0.04, 'mono', 6.3E8, 120,
                        None, 'oa', 'Citric acid', 'mono']
    df04 = out[out['sample_time'] < 213]
    fig = px.line(df04, x='sample_time', category_orders={'media': ['Citric acid', 'Citric acid + Thiamine']} hover_data=[
                'sample_time'], y='CFUs', log_y=True, title='ct_04', facet_col='media', color='species', line_dash='mono')

    df15 = out[(out['sample_time'] > 213) & (out['sample_time'] < 381)]
    fig = px.line(df15, x='sample_time', hover_data=[
                'sample_time'], y='CFUs', log_y=True, title='ct_04', facet_col='media', color='species', line_dash='mono')

    df085 = out[out['sample_time'] > 381]
    fig = px.line(df085, x='sample_time', hover_data=[
                'sample_time'], y='CFUs', log_y=True, title='ct_04', facet_col='media', color='species', line_dash='mono')

    figs.append(fig)
    return figs

dilution_rates = {0.04: [45.5, 212.97], 0.15: [
    212.97, 381.466], 0.085: [381.466, 619], 'all': [45.4, 619]}
reactors = {'M0': ['ct'], 'M2': ['ct'], 'M4': ['ct', 'oa'], 'M5': ['ct', 'oa']}

def dump_dfs():
    od = chibio_parser('citrate_thiamine_merged', down_sample=True)[0]
    cfus = cfu_parser('citrate_thiamine_merged')[0]
    for r, species in reactors.items():
        d = join('dfs', r)
        if not exists(d):
            mkdir(d)
        for dr, period in dilution_rates.items():
            d = join('dfs', r, str(dr))
            if not exists(d):
                mkdir(d)
            mask = (od['reactor'] == r) & (od['sensor'] == 'od_measured') \
                & (od['exp_time'] >= period[0]) & (od['exp_time'] <= period[1])
            od[mask].to_csv(join('dfs', r, str(dr), 'od.csv'), index=False)
            for s in species:
                mask = (cfus['reactor'] == r) & (cfus['species'] == s) \
                    & (cfus['sample_time'] >= period[0]) & (cfus['sample_time'] <= period[1])
                cfus[mask].to_csv(
                    join('dfs', r, str(dr), s+'.csv'), index=False)


def plot():
    base_dir = join('/', 'home', 'eric', 'ChiBioFlow',
                    'data', 'citrate_thiamine_merged')
    colors = ['#8872cd', '#e27e50']
    names = {'ct': 'Ct',
             'oa': 'Oa',
             'OD': 'OD'}
    title = ['Ct, C + T', 'Ct, C', 'Ct + Oa, C', 'Ct + Oa, C + T']
    figs = []
    for dr, period in dilution_rates.items():
        y_od = []
        y_cfus = []
        fig = make_subplots(rows=2, cols=len(reactors.keys()))
        for j, (r, species) in enumerate(reactors.items()):
            dr = str(dr)
            od = pd.read_csv(join(base_dir, 'dfs', r, dr, 'od.csv'))
            fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                                     marker=dict(color='#454242'), name='OD', mode="lines", yaxis='y1'), row=1, col=1+j)
            y_od.append(fig['data'][-1]['yaxis'])
            for i, s in enumerate(species):
                cfus = pd.read_csv(join(base_dir, 'dfs', r, dr, s + '.csv'))
                fig.add_trace(go.Scatter(x=cfus['sample_time'], error_y=dict(type='data', array=cfus['stdev'], visible=True),
                                         y=cfus['average'], marker=dict(color=colors[i]), name=s), row=2, col=1+j)
                y_cfus.append(fig['data'][-1]['yaxis'])
        for axis in fig['layout']:
            if axis in ['yaxis', 'yaxis2', 'yaxis3', 'yaxis4']:
                fig['layout'][axis]['range'] = [0, 0.6]
            if axis in ['yaxis5', 'yaxis6', 'yaxis7', 'yaxis8']:
                fig['layout'][axis]['type'] = 'log'
                fig['layout'][axis]['range'] = [6, 10]
        fig.update_xaxes(dict(range=[period[0], period[1]]))
        legend = []
        for d in fig['data']:
            d['name'] = names[d['name']]
            if d['name'] not in legend:
                d['showlegend'] = True
            else:
                d['showlegend'] = False
            legend.append(d['name'])
        ts = set(cfus['sample_time'])
        for t in ts:
            fig.add_vline(x=t, opacity=0.3)
        fig.show()
        figs.append(fig)
    return figs
