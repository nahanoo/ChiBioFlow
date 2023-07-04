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
    out = pd.DataFrame(columns=['D', 'cond', 'CFUs', 'sample_time'])
    df, o = cfu_parser('citrate_thiamine_merged')
    # df.astype({'sample_time':'float'})
    df = df[df['sample_time'] < 213]
    conditions = {'M0': 'mono_thiamine', 'M2': 'mono',
                  'M4': 'comm', 'M5': 'comm_thiamine'}
    for r, s, c, t in zip(df['reactor'], df['species'], df['average'], df['sample_time']):
        if (r in conditions.keys()) & (s == 'ct'):
            out.loc[len(out)] = [0.04, conditions[r], c, t]

    fig = px.box(out, x='cond', hover_data=['sample_time'], y='CFUs', log_y=True, category_orders={
        'cond': ['mono', 'comm', 'mono_thiamine', 'comm_thiamine']}, points='all', height=250, width=400)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.show()

    out = pd.DataFrame(columns=['D', 'cond', 'CFUs'])
    df, o = cfu_parser('citrate_thiamine_merged')
    # df.astype({'sample_time':'float'})
    df = df[df['sample_time'] > 213]
    conditions = {'M0': 'mono_thiamine', 'M2': 'mono',
                  'M4': 'comm', 'M5': 'comm_thiamine'}
    for r, s, c in zip(df['reactor'], df['species'], df['average']):
        if (r in conditions.keys()) & (s == 'ct'):
            out.loc[len(out)] = [0.04, conditions[r], c]

    fig = px.box(out, x='cond', y='CFUs', log_y=True,  category_orders={
        'cond': ['mono', 'comm', 'mono_thiamine', 'comm_thiamine']}, points='all', height=250, width=400)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.show()


titles = ['Oa', 'Ct + Oa', 'Ct + Oa + Thiamine',
          'Ct + Thiamine', 'Oa + Thiamine', 'Ct']


dilution_rates = {0.04: [45.5, 212.97], 0.15: [
    212.97, 381.466], 0.085: [381.466, 481.89], 'all': [45.4, 481.89]}
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
    colors = ['#8872cd', '#e27e50']
    names = {'ct':'Ct',
           'oa':'Oa',
           'OD':'OD'}
    for dr, period in dilution_rates.items():
        y_od = []
        y_cfus = []
        fig = make_subplots(rows=2, cols=len(reactors.keys()))
        for j, (r, species) in enumerate(reactors.items()):
            dr = str(dr)
            od = pd.read_csv(join('dfs', r, dr, 'od.csv'))
            fig.add_trace(go.Scatter(x=od['exp_time'], y=od['measurement'],
                                     marker=dict(color='#454242'), name='OD', mode="lines", yaxis='y1'), row=1, col=1+j)
            y_od.append(fig['data'][-1]['yaxis'])
            for i, s in enumerate(species):
                cfus = pd.read_csv(join('dfs', r, dr, s + '.csv'))
                fig.add_trace(go.Scatter(x=cfus['sample_time'],
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
        return fig
