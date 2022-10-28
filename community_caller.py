import plotly.express as px
import plotly.graph_objects as go
from scipy.integrate import odeint
import numpy as np
from model_community import Chain
from plotting import plot_community_model, plot_dilution_factors, plot_carrying_capacity
from math import exp


def update_labels(fig, x_label, y_label, title, size=[250, 300]):
    fig.for_each_xaxis(lambda axis: axis.title.update(text=x_label))
    fig.for_each_xaxis(lambda axis: axis.title.update(text=x_label))
    #fig.for_each_yaxis(lambda axis: axis.title.update(text=y_label))
    fig.update_layout(yaxis_title=y_label, title=title,
                      height=size[0], width=size[1], paper_bgcolor='#999999', plot_bgcolor='#999999', font_color='black')
    fig.update_xaxes(showline=True, linewidth=0.5,
                     linecolor='black', gridwidth=0.5, gridcolor='black', zerolinecolor='black')
    fig.update_yaxes(showline=True, linewidth=0.5,
                     linecolor='black', gridwidth=0.5, gridcolor='black', zerolinecolor='black', zerolinewidth=0.5)
    for anno in fig['layout']['annotations']:
        anno['text'] = ''
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
    xs = range(t)
    ys = [s[0] for s in odeint(logistic_model, 0.08, xs, args=(k, K))]
    fig = px.line(x=xs, y=ys)
    fig = update_labels(fig, 'Time in hours', 'OD', "Logistic model")
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


def plot_chain(t, D, n_reactors):
    chain = get_chain(t, D, n_reactors)
    fig = plot_community_model(chain)
    fig = update_labels(fig, 'Time in hours', 'OD', "$D < \mu$")
    return fig


def plot_capacity(t, n_reactors):
    D = 0.32899
    chain = get_chain(t, D, n_reactors, capacity=True)
    fig = plot_carrying_capacity(chain)
    fig = update_labels(fig, 'Time in hours', 'K', "Carrying capacity",size=[250,700])
    return fig


def plot_lim_growth_rates(t):
    xs = np.arange(0, t, 0.1)
    ys = ys = [s[0]
               for s in odeint(logistic_model, 0.08, xs, args=(0.303129065, 1.5))]
    r_logs = [logistic_model(y, None, 0.303129065, 1.5)/y for y in ys]
    rs = [exp(r_log) - 1 for r_log in r_logs]
    fig = px.line(x=ys, y=rs)
    fig.update_layout(
        xaxis=dict(
            tickmode='linear',
            dtick=0.2
        ))
    fig = update_labels(fig, 'N', '$\mu$', '$\mu = f(N)$', size=[250, 300])
    return fig


def plot_growth_rates(t):
    interval = 0.05
    period = int(1 / interval)
    xs = np.arange(0, t, interval)
    ys = [s[0]
          for s in odeint(logistic_model, 0.00001, xs, args=(0.303129065, 1.5))]
    us = []
    for i, x in enumerate(xs):
        try:
            us.append((ys[i + period] - ys[i]) / ys[i])
        except IndexError:
            break

    fig = px.line(x=ys[:-period], y=us)
    fig.update_layout(
        xaxis=dict(
            tickmode='linear',
            dtick=0.2
        ))
    fig = update_labels(fig, 'N', '$\mu$', '$\mu = f(N)$', size=[250, 300])
    return fig
