import plotly.express as px
import plotly.graph_objects as go
from scipy.integrate import odeint
import numpy as np
from model_community import Chain
from plotting import plot_composition, plot_community_model, plot_species, plot_od, plot_dilution_factors, plot_carrying_capacity, plot_species_model, style_plot, plot_species_model_composition, plot_od
from math import exp
from model_species import Chain as Species_chain
import pandas as pd
from species_main import random_compositions
from fitting import fit_comm
from samples import Samples


def update_labels(fig, x_label, y_label, title, size=[250, 300], width=6, style_colors=True, font_size=18, reactor_tiles=False):
    fig.for_each_xaxis(lambda axis: axis.title.update(text=None))
    fig.for_each_yaxis(lambda axis: axis.title.update(text=None))
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
    fig.update_layout(font={'size': font_size})
    fig.update_traces(line={'width': width})
    if style_colors:
        fig.update_traces(line={'color': line_color})
    fig['layout']['legend']['title']['text'] = ''
    if reactor_tiles:
        titles = ['Reactor '+str(n) for n in range(1, 5)]
        for counter, annotation in enumerate(fig['layout']['annotations']):
            annotation['text'] = titles[counter]

    return fig


def exp_model(N, t, k):
    return k * N


def logistic_model(N, t, k, K):
    return k * N * (1 - N/K)


def plot_dilution(t, n_reactors):
    D = 0.3
    chain = get_chain(t, D, n_reactors, transfer_rate=1)
    fig = plot_dilution_factors(chain)[1]
    fig = update_labels(fig, 'Time [h]', '$\large{D_{Effective}}$', "", size=[
                        200, 450], width=3, font_size=16, reactor_tiles=True)
    for i, d in enumerate(fig.data):
        Fd = str(chain.chain[i].dilution_factors[-1])[:4]
        fig.add_trace(go.Scatter(x=[20], y=[Fd],
                                 mode='markers+text',
                                 text=Fd,
                                 textposition='middle right',
                                 marker=dict(color=d.line.color, size=12),
                                 legendgroup=d.name,
                                 showlegend=False), row=1, col=i + 1)
    fig.update_layout(
        margin=dict(l=100, r=10, t=50, b=10),
    )
    # fig.update_yaxes(dtick=0.25)
    # fig.update_xaxes(dtick=5)
    return fig


def plot_exp_chain():
    fig = plot_od('c1')[1]
    fig = update_labels(fig, 'Time [h]', 'OD', '', size=[
                        250, 750], width=2, reactor_tiles=True)
    fig.update_layout(showlegend=False)
    return fig


def plot_logistic_growth(t, k, K):
    xs = np.arange(0, t, 0.1)
    ys = [s[0] for s in odeint(logistic_model, 0.01, xs, args=(k, K))]
    fig = px.line(x=xs, y=ys)
    fig = update_labels(
        fig, 'Time [h]', 'OD', "$\large{r_{max} = 0.28, K = 1.5}$", size=[250, 250])
    return fig


def get_chain(t, D, n_reactors, transfer_rate=1):
    # Steady state D = 0.32899
    chain = Chain(n_reactors)
    chain.transfer_rate = transfer_rate
    chain.chain[0].N = 0.1
    chain.dilution_rate = D
    chain.experiment(t)
    return chain


def gradients():
    Ds = [0.1, 0.2, 0.3, 0.4]
    Df = pd.DataFrame(columns=['reactor', 'Df', 'D'])
    Ks = pd.DataFrame(columns=['reactor', 'Ks', 'D'])
    Ns = pd.DataFrame(columns=['reactor', 'Ns', 'D'])
    i = 0
    figs = []
    for D in Ds:
        chain = get_chain(60, D, 4)
        for c in chain.chain:
            Df.loc[i] = [c.name + 1, c.dilution_factors[-1], D]
            Ks.loc[i] = [c.name + 1, c.Ks[-1], D]
            Ns.loc[i] = [c.name + 1, c.N, D]
            i += 1
    fig = px.line(Df, x='reactor', y='Df', color='D', line_dash='D')
    fig = update_labels(
        fig, 'Reactor', '$\large{D_{Effective}}$', '$\large{D_{Effective}}\space gradient$')
    fig.update_xaxes(type='category')
    figs.append(fig)
    fig = px.line(Ks, x='reactor', y='Ks', color='D', line_dash='D')
    fig = update_labels(fig, 'Reactor', 'K', 'K gradient')
    fig.update_xaxes(type='category')
    figs.append(fig)
    fig = px.line(Ns, x='reactor', y='Ns', color='D', line_dash='D')
    fig = update_labels(fig, 'Reactor', 'N', 'N gradient')
    fig.update_xaxes(type='category')
    figs.append(fig)
    return figs


def get_species_chain(t, n_reactors, D=0.31395, capacity=False, transfer_rate=0):
    # Steady state D = 0.32899
    chain = Species_chain(n_reactors)
    chain.transfer_rate = transfer_rate
    chain.dilution_rate = D

    for specie in chain.chain[0].species.values():
        specie.N = 0.025

    for c in chain.chain[1:]:
        for specie in c.species.values():
            specie.N = 0
    chain.experiment(t)
    return chain


def plot_species_chain(t, D, n_reactors):
    chain = get_species_chain(t, n_reactors, transfer_rate=1, D=D)
    fig = plot_species_model(chain)
    fig = update_labels(fig, 'Time [h]', 'OD',
                        'Species model', size=[250, 500], style_colors=False)
    return fig


def plot_chain(t, D, n_reactors):
    chain = get_chain(t, D, n_reactors, transfer_rate=2)
    fig = plot_community_model(chain)
    fig = update_labels(fig, 'Time in hours', 'OD',
                        "$D < \mu$", reactor_tiles=True)
    return fig


def plot_capacity(t, n_reactors):
    D = 0.32899
    chain = get_chain(t, D, n_reactors)
    fig = plot_carrying_capacity(chain)
    fig = update_labels(fig, 'Time [h]', 'K',
                        "D = 0.32", size=[250, 700])
    return fig


def Nt(t, N, k, K):
    xs = [0, t]
    return odeint(logistic_model, N, xs, args=(k, K))[-1][0]


def specific_growth_rate(t, k, K, N):
    f = 1/Nt(t, N, k, K)
    dNdt = k * Nt(t, N, k, K) * (1 - Nt(t, N, k, K)/K)
    return (exp(f * dNdt) - 1)


def growth_rate_community(t, K=1.5):
    xs = np.arange(0, t, 0.5)
    Ns = [N[0]
          for N in odeint(logistic_model, 0.1, xs, args=(0.28, K))]
    rs = [specific_growth_rate(t, 0.28, K, 0.1) for t in xs]
    fig = px.line(x=Ns, y=rs)
    fig = update_labels(
        fig, 'OD', '$\large{\mu\space [h^{-1}]}$', 'Community')
    fig.update_layout(
        margin=dict(l=100, r=10, t=50, b=10),
    )
    fig.update_xaxes(dtick=0.2, tick0=0.1)
    return fig


def growth_rate_species(t):
    chain = get_species_chain(t, 1, transfer_rate=0)
    dfs = []
    xs = chain.xs
    e = 'species_growth_rates'
    cfus, fig = plot_species(e)
    for c in chain.chain:
        rates = {name: None for name in c.species.keys()}
        for specie, params in c.species.items():
            s_df = cfus[cfus['species'] == specie]
            t0_df = s_df[s_df['sample_time'] == 0.0]
            N0 = t0_df.loc[t0_df.index[0]]['average']
            t70_df = s_df[s_df['sample_time'] == 25.0]
            N70 = t70_df.loc[t70_df.index[0]]['average']
            Ns = params.ys
            rates[specie] = [specific_growth_rate(
                t, params.r, N70, N0) for t in xs]
            df = pd.DataFrame(
                {'rate': rates[specie], 'specie': specie, 'N': c.total})
            dfs.append(df)

    df = pd.concat(dfs)
    fig = px.line(df, x='N', y='rate', color='specie')
    fig = style_plot(None, fig, 'cfus')
    fig = update_labels(fig, 'OD', '$\mu\space [h^{-1}]$',
                        'Species', size=[250, 450], style_colors=False, font_size=18)
    fig.update_layout(
        margin=dict(l=100, r=10, t=50, b=10),
    )
    return fig


def plot_species_chain(t, D, n_reactors):
    chain = get_species_chain(t, n_reactors, D=D, transfer_rate=1)
    fig = plot_species_model(chain)
    fig = style_plot(None, fig, 'cfus')
    fig = update_labels(fig, 'Time [h]', 'OD 600nm',
                        str(D), size=[250, 750], style_colors=False)

    return fig


def random_compositions(t, compositions):
    dfs = []
    compositions = [np.random.dirichlet(np.ones(4), size=1)[
        0] for i in range(compositions)]
    for counter, composition in enumerate(compositions):
        chain = Species_chain(4)
        chain.dilution_rate = 0.15
        chain.transfer_rate = 2
        for j, specie in enumerate(chain.chain[0].species.values()):
            specie.N = np.random.choice(
                np.arange(0.025, 0.5, 0.1)) * compositions[counter][j]

        for c in chain.chain[1:]:
            for specie in c.species.values():
                specie.N = 0

        chain.experiment(t)

        for c in chain.chain:
            for name, specie in c.species.items():
                df = pd.DataFrame({'x': chain.xs, 'N': specie.ys,
                                   'species': name, 'reactor': c.name, 'composition': counter})
                dfs.append(df)

    out = pd.concat(dfs)
    fig = px.line(out, x='x', y='N', facet_col='reactor',
                  color='species', line_group='composition', facet_col_spacing=0.05)

    fig = style_plot(None, fig, 'cfus')
    fig = update_labels(fig, 'Time [h]', 'OD',
                        'Species model D = 0.15', size=[250, 750], width=2)

    return fig


def species_composition(t):
    chain = Species_chain(4)
    chain.transfer_rate = 1
    chain.dilution_rate = 0.15

    chain.chain[0].species['at'].N = 0.037
    chain.chain[0].species['ct'].N = 0.0066
    chain.chain[0].species['ms'].N = 0.043
    chain.chain[0].species['oa'].N = 0.013
    for c in chain.chain[1:]:
        for specie in c.species.values():
            specie.N = 0
    chain.experiment(72)
    fig, df = plot_species_model_composition(chain)
    fig = style_plot(None, fig, 'cfus')
    fig = update_labels(
        fig, 'Time [h]', 'Composition', None, size=[200, 750], font_size=16, style_colors=False, width=3)
    return fig


def species_composition_D():
    t = 72
    ds = np.arange(0.05, 0.45, 0.05)
    dfs = []
    for d in ds:
        chain = get_species_chain(t, 4, transfer_rate=1, D=d)
        df = plot_species_model_composition(chain)[1]
        df.insert(len(df.columns), 'Dilution rate', d)
        dfs.append(df)
    out = pd.concat(dfs)
    out = out[out['x'] == t]
    fig = px.line(out, x='Dilution rate', y='ratio', color='species',
                  facet_col='reactor', facet_col_spacing=0.05)
    fig = style_plot(None, fig, 'cfus')
    fig = update_labels(fig, '$D\space h^{-1}$', 'Composition',
                        None, size=[250, 750], style_colors=False, width=3, reactor_tiles=True)
    return fig


def log_community():
    fig = fit_comm()
    fig = update_labels(fig, 'Time [h]', 'OD',
                        "$\large{r_{max} = 0.28, K = 1.5}$", width=3)
    colors = ['#372f55', '#7858a6']
    for i, line in enumerate(fig['data']):
        line['line']['color'] = colors[i]
    fig['layout']['legend']['title']['text'] = ''
    return fig


def plot_experiment_model():
    e = 'c1'
    chain = get_chain(52, 0.3, 4, transfer_rate=2)
    df, fig = plot_od(e, model=True, chain=chain)
    fig = update_labels(fig, 'Time [h]', 'OD',
                        'D = 0.3', size=[250, 750], width=2, reactor_tiles=True)
    for data in fig['data']:
        if data['name'] == 'Model':
            data['name'] = 'model'
            data['line']['color'] = '#372f55'
        if data['name'] == 'OD measured':
            data['name'] = 'OD'
        if (data['name'] == 'OD') | (data['name'] == ''):
            data['line']['color'] = '#7858a6'
    return fig


def plot_species_rates():
    e = 'species_growth_rates'
    df, fig = plot_species(e)
    fig = update_labels(fig, 'Time [h]', 'CFUs/mL', 'Species growth curves',
                        style_colors=False, size=[250, 400], width=3)
    return fig


def species_fit():
    chain = Species_chain(1)
    e = 'species_growth_rates'
    df, fig = plot_species(e)
    fig.update_traces(line={'dash': 'dot'})
    t = 70
    s = Samples()
    c = chain.chain[0]
    out = pd.DataFrame(columns=['at', 'ct', 'ms', 'oa', 'time'])
    xs = np.arange(0, t + 0.5, 0.5)
    out['time'] = xs
    for name, specie in c.species.items():
        s_df = df[df['species'] == name]
        t0_df = s_df[s_df['sample_time'] == 0.0]
        N0 = t0_df.loc[t0_df.index[0]]['average']
        if name == 'ms':
            t70_df = s_df[s_df['sample_time'] == 20.0]
            N70 = t70_df.loc[t70_df.index[0]]['average']
        else:
            t70_df = s_df[s_df['sample_time'] == 70.0]
            N70 = t70_df.loc[t70_df.index[0]]['average']

        ys = [s[0]
              for s in odeint(logistic_model, N0, xs, args=(specie.r, N70))]
        out[name] = ys
    colors = {'at': '#2c8c5a',
              'ct': '#8872cd',
              'oa': '#e27e50',
              'ms': '#e5b008'}
    for counter, f in enumerate(fig['data']):
        name = list(c.species.keys())[counter]
        fig.add_trace(go.Scatter(x=out['time'], y=out[name], showlegend=False, opacity=0.8, line={
                      'color': colors[name]}))
    fig.add_trace(go.Scatter(x=[0], name='model', y=None,
                  line={'color': 'black'}))
    fig['layout']['legend']['title']['text'] = ''
    fig = update_labels(fig, 'Time [h]', 'CFUs/mL', '',
                        style_colors=False, size=[250, 400], width=3, font_size=16)
    return df, fig


def plot_exp_species_comp():
    e = 'd_0.15'
    df, fig = plot_composition(e)
    fig = update_labels(fig, 'Time [h]', 'Composition', None, size=[
                        200, 750], style_colors=False, font_size=16, width=3, reactor_tiles=True)
    fig.update_traces(line={'dash': 'dot'})
    return fig


def plot_cfus_d_028():
    e = 'd_0.15'
    df, fig = plot_species(e)
    fig = update_labels(fig, 'Time [h]', 'CFUS/mL', None, size=[250, 750],
                        style_colors=False, width=3, font_size=16, reactor_tiles=True)
    return fig
