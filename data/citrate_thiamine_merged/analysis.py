import plotly.express as px
import pandas as pd
from glob import glob
from os.path import join, split
from chibio_parser import *
from plotting import *


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
    out = pd.DataFrame(columns=['D', 'cond', 'CFUs','sample_time'])
    df, o = cfu_parser('citrate_thiamine_merged')
    # df.astype({'sample_time':'float'})
    df = df[df['sample_time'] < 213]
    conditions = {'M0': 'mono_thiamine', 'M2': 'mono',
                  'M4': 'comm', 'M5': 'comm_thiamine'}
    for r, s, c,t in zip(df['reactor'], df['species'], df['average'],df['sample_time']):
        if (r in conditions.keys()) & (s == 'ct'):
            out.loc[len(out)] = [0.04, conditions[r], c,t]

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


def plot_od():
    df = chibio_parser('citrate_thiamine_merged', down_sample=True)[0]
    df = df[(df['measurement'] <= 0.7) & (df['sensor'] == 'od_measured')]
    df = df.sort_values(by='reactor')

    fig = plot_od('citrate_thiamine_merged', df=df)[1]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]


def species():
    fig = plot_species('citrate_thiamine_merged')[1]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.show()

    fig = plot_composition('citrate_thiamine_merged')[1]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.show()
