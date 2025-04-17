import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from glob import glob
import pandas as pd
from plate_reader import *
from ct_chibio_batch import *
from oa_chibio_batch import *


colors = {"ct": "#7570B3", "oa": "#D95F02", "total": "#565e91"}

w = 400
h = 250


def style_plot(fig, w, h, marker_size=3):
    """Style function for figures setting fot size and true black color."""
    for d in fig["data"]:
        d["marker"]["size"] = marker_size
        d["line"]["width"] = marker_size
    # Font size
    j = 10
    fig.update_layout(font={"size": j, "color": "black"})
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = j
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = j
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = j
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"
    fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=j, color="black")))
    fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(size=j, color="black")))
    fig.update_layout(
        margin=dict(l=60, r=20, t=30, b=20),
    )
    fig.update_yaxes(title_standoff=10)
    fig.update_xaxes(title_standoff=10)
    fig.update_layout(width=w, height=h)
    return fig


def plot_gradient(fit=False):
    figs = plot_plate_reader(fit=fit)
    for fig in figs:
        fig = style_plot(fig, w * 2, h * 2, marker_size=2)
        fig.show()


def plot_growth_rates_plate_reader():
    figs = plate_reader_growth_rates()
    for fig in figs:
        fig = style_plot(fig, w * 2, h * 1.5, marker_size=2)
        fig.show()


def plot_max_OD_plate_reader():
    figs = plate_reader_max_OD()
    for fig in figs:
        fig = style_plot(fig, w * 2, h * 1.5, marker_size=2)
        fig.show()


def plot_chibio_batch(plot_log=True):
    fig = make_subplots(
        rows=1,
        cols=2,
        column_titles=["Comamonas testosteroni", "Ochrobactrum anthropi"],
        shared_yaxes=True,
        shared_xaxes=True,
    )
    fig_ct = plot_ct_chibio_batch(plot_log=plot_log)
    fig_oa = plot_oa_chibio_batch(plot_log=plot_log)
    for trace in fig_ct.data:
        fig.add_trace(trace, row=1, col=1)

    for trace in fig_oa.data:
        fig.add_trace(trace, row=1, col=2)
    fig = style_plot(fig, w * 3, h * 1.5, marker_size=2)
    fig.show()
