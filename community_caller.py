import plotly.express as px
from scipy.integrate import odeint
import numpy as np
from model_community import Chain
from plotting import plot_community_model, plot_dilution_factors, plot_carrying_capacity, plot_species_model,style_plot
from math import exp
from model_species import Chain as Species_chain
import pandas as pd


def update_labels(fig, x_label, y_label, title, size=[250, 300]):
    fig.for_each_xaxis(lambda axis: axis.title.update(text=x_label))
    fig.for_each_xaxis(lambda axis: axis.title.update(text=x_label))
    # fig.for_each_yaxis(lambda axis: axis.title.update(text=y_label))
    lw = 1
    color = '#999999'
    fig.update_layout(yaxis_title=y_label, title=title,
                      height=size[0], width=size[1], paper_bgcolor=color, plot_bgcolor=color, font_color='black')
    fig.update_xaxes(showline=True, linewidth=lw,
                     linecolor='black', gridwidth=lw, gridcolor='black', zerolinecolor='black')
    fig.update_yaxes(showline=True, linewidth=lw,
                     linecolor='black', gridwidth=lw, gridcolor='black', zerolinecolor='black', zerolinewidth=lw)
    for anno in fig['layout']['annotations']:
        anno['text'] = ''
    # fig.update_layout(title=None)
    fig.update_layout(
        margin=dict(l=10, r=10, t=50, b=10),
    )
    fig.update_layout(font={'size': 18})
    fig.update_traces(line={'width': 6})
    return fig


def exp_model(N, t, k):
    return k * N


def logistic_model(N, t, k, K):
    return k * N * (1 - N/K)


def plot_dilution(t, n_reactors):
    D = 0.32899
    chain = get_chain(t, D, n_reactors)
    fig = plot_dilution_factors(chain)
    fig = update_labels(fig, 'Time in hours',
                        'Dilution factor', 'Dilution factors')
    return fig


def plot_logistic_growth(t, k, K):
    xs = np.arange(0, t, 0.1)
    ys = [s[0] for s in odeint(logistic_model, 0.01, xs, args=(k, K))]
    fig = px.line(x=xs, y=ys)
    fig = update_labels(fig, 'Time [h]', 'OD 600nm', "Logistic model")
    return fig


def get_chain(t, D, n_reactors, capacity=False):
    # Steady state D = 0.32899
    chain = Chain(n_reactors)
    chain.transfer_rate = 1
    chain.chain[0].N = 0.08
    chain.dilution_rate = D
    if capacity:
        chain.reduce_carrying_capacity = True
    chain.experiment(t)
    return chain


def get_species_chain(t, D, n_reactors, capacity=False):
    # Steady state D = 0.32899
    chain = Species_chain(n_reactors)
    chain.transfer_rate = 0
    for specie in chain.chain[0].species.values():
        specie.N = 0.025

    for c in chain.chain[1:]:
        for specie in c.species.values():
            specie.N = 0
    chain.experiment(t)
    return chain


def plot_species_chain(t, D, n_reactors):
    chain = get_species_chain(t, D, n_reactors)
    fig = plot_species_model(chain)
    fig = update_labels(fig, 'Time [h]', 'OD 600nm',
                        'Species model', size=[250, 500])
    return fig


def plot_chain(t, D, n_reactors):
    chain = get_chain(t, D, n_reactors)
    chain.chain[0].N = 0.1
    fig = plot_community_model(chain)
    fig = update_labels(fig, 'Time in hours', 'OD', "$D < \mu$")
    return fig


def plot_capacity(t, n_reactors):
    D = 0.32899
    chain = get_chain(t, D, n_reactors, capacity=True)
    fig = plot_carrying_capacity(chain)
    fig = update_labels(fig, 'Time in hours', 'K',
                        "Carrying capacity", size=[250, 700])
    return fig


def Nt(t, N,k,K):
    xs = [0, t]
    return odeint(logistic_model, N, xs, args=(k, K))[-1][0]


def specific_growth_rate(t, k, K, N):
    f = 1/Nt(t, N,k,K)
    dNdt = k * Nt(t, N,k,K) * (1 - Nt(t, N,k,K)/K)
    return (exp(f * dNdt) - 1)


def growth_rate_community(t):
    xs = np.arange(0, t, 0.1)
    Ns = [N[0]
          for N in odeint(logistic_model, 0.1, xs, args=(0.303129065, 1.5))]
    rs = [specific_growth_rate(t, 0.303129065, 1.5, 0.1) for t in xs]
    fig = px.line(x=Ns, y=rs)
    fig = update_labels(fig, 'N', '$\mu$', 'Specific growth rate')
    return fig

def growth_rate_species(t):
    chain = get_species_chain(60, None, 1)
    c = chain.chain[0]
    rates = {name: None for name in c.species.keys()}
    df = pd.DataFrame()
    xs = np.arange(0, t, 0.5)
    dfs = []
    for specie, params in c.species.items():
        Ns = [N[0]
            for N in odeint(logistic_model, 0.025, xs, args=(params.r, params.K * 1.5))]
        rates[specie] = [specific_growth_rate(
            t, params.r, params.K * 1.5, 0.025) for t in xs]
        df = pd.DataFrame({'rate':rates[specie],'specie':specie,'N':Ns})
        dfs.append(df)

    df = pd.concat(dfs)
    fig = px.line(df,x='N',y='rate',color='specie')
    fig = style_plot(None,fig,'cfus')
    fig = update_labels(fig,'OD 600nm','$\mu$','Species specific growth rates',size=[250, 500])
    return fig