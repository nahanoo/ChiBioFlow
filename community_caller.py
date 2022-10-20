import plotly.express as px
from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
from model_community import Chain
from plotting import plot_community_model

def update_labels(fig, x_label, y_label,title):
    fig.for_each_xaxis(lambda axis: axis.title.update(text=x_label))
    fig.for_each_yaxis(lambda axis: axis.title.update(text=y_label))
    fig.update_layout(title=title)
    return fig


def exp_model(N, t, k):
    return k * N


def plot_exp_growth(t, k):
    xs = range(t)
    ys = [s[0] for s in odeint(exp_model, 0.08, xs, args=(k,))]
    fig = px.line(x=xs, y=ys)
    fig = update_labels(fig,'Time in hours','OD',"Exponential growth")
    fig.update_layout(height=250)
    fig.update_layout(width=300)
    fig.show()

def plot_reactor_increase(t):
    chain = Chain(1)
    chain.model = 'exponential'
    chain.transfer_rate = 1
    chain.chain[0].N = 0.08
    chain.dilution_rate = 0.15
    chain.experiment(t)
    fig = plot_community_model(chain)
    fig = update_labels(fig,'Time in hours','OD',"D < u")
    fig.show()


def plot_reactor_steady(t):
    chain = Chain(1)
    chain.model = 'exponential'
    chain.transfer_rate = 1
    chain.chain[0].N = 0.08
    chain.dilution_rate = 0.35408923072260035
    chain.experiment(t)
    fig = plot_community_model(chain)
    fig.show()


