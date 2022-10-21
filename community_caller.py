import plotly.express as px
from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
from model_community import Chain
from plotting import plot_community_model


def update_labels(fig, x_label, y_label, title, size=[250, 300]):
    fig.for_each_xaxis(lambda axis: axis.title.update(text=x_label))
    fig.for_each_yaxis(lambda axis: axis.title.update(text=y_label))
    fig.update_layout(title=title)
    fig.update_layout(height=size[0])
    fig.update_layout(width=size[1])
    for anno in fig['layout']['annotations']:
        anno['text'] = ''
    return fig


def exp_model(N, t, k):
    return k * N


def logistic_model(N, t, k, K):
    return k * N * (1 - N/K)


def plot_logistic_growth(t, k, K):
    xs = range(t)
    ys = [s[0] for s in odeint(logistic_model, 0.08, xs, args=(k, K))]
    fig = px.line(x=xs, y=ys)
    fig = update_labels(fig, 'Time in hours', 'OD', "Logistic model")
    return fig


def plot_chain(t, D, n_reactors):
    # Steady state D = 0.32899
    chain = Chain(n_reactors)
    chain.transfer_rate = 1
    chain.chain[0].N = 0.08
    chain.dilution_rate = D
    chain.experiment(t)
    fig = plot_community_model(chain)
    fig = update_labels(fig, 'Time in hours', 'OD', "$D < \mu$")
    return fig
