from genericpath import exists
import plotly.express as px
import argparse
from os.path import join
import pandas as pd
from glob import glob
from chibio_parser import cfu_parser
from chibio_parser import chibio_parser
import numpy as np


def parse_args():
    """Parsing variables for plotting."""
    parser = argparse.ArgumentParser(
        description="Plotting library for ChiBio.")
    parser.add_argument("experiment", help="name of the experiment directory")
    parser.add_argument("mode", help="use 'chibio' for plotting a column of the data files; \
        use 'strain' to plot strain composition.")
    parser.add_argument(
        "--column", help="column name to plot from ChiBio csv data")

    return parser.parse_args()


def main():
    global e
    global c
    global chain

    args = parse_args()
    e = args.experiment
    c = args.column
    mode = args.mode
    with open(join('data', e, 'order.txt'), 'r') as f:
        chain = f.read().rstrip().split(',')

    if mode == 'chibio':
        df, fig = plot_chibio()
        fig = style_plot(fig, 'od_measured')
        fig.show()
        return df, fig

    if mode == 'species':
        df = cfu_parser(e, chain)
        fig = plot_species(df)
        fig = style_plot(fig, 'cfus')
        fig.show()
        return df, fig

    if mode == 'composition':
        df = cfu_parser(e, chain)
        fig = plot_composition(df)
        fig = style_plot(fig, 'composition')
        fig.show()
        return df, fig

    if mode == 'total':
        df = cfu_parser(e, chain)
        fig = plot_total(df)
        fig = style_plot(fig, 'cfus')
        fig.show()
        return df, fig

    return df, fig


def plot_chibio():
    """Creates lineplot for parsed column e.g. od_measured.
    Plots every reactor as subplot.
    """
    
    df = chibio_parser(e,chain,c)
    fig = px.line(df, x="exp_time", y=c, facet_col="reactor",
                  facet_col_wrap=2, hover_data=['exp_time'])

    return df, fig


def plot_total(df):
    """Plot sum of CFUs of all species"""
    df = df[['reactor','sample_time','total']].drop_duplicates()
    fig = px.line(df, x="sample_time", y='total', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, log_y=True)
    fig.update_layout(font={'size': 20},
                      xaxis_title='Time in hours',
                      yaxis_title='CFUs/mL')
    return fig


def plot_species(df):
    """Plots CFUs based on parsed xlsx sheet"""
    df['average'] = df['average'].replace(0, np.nan)
    fig = px.line(df, x="sample_time", y='average', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, error_y='stdev', color='species', log_y=True)
    return fig


def plot_composition(df):
    """Plots community composition in percent"""
    fig = px.line(df, x="sample_time", y='composition', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, color='species', log_y=False)
    return fig


def style_plot(fig, style, fontsize=20):
    """Updated labels and titles."""
    def species_colors(fig):
        # Species color code
        colors = {'at': '#2c8c5a',
                  'ct': '#8872cd',
                  'oa': '#e27e50',
                  'ms': '#e5b008'}
        for data in fig['data']:
            try:
                data['line']['color'] = colors[data['name']]
            except KeyError:
                pass
        return fig

    def species_names(fig):
        # Species names
        names = {'at': '<i>A. tumefaciens</i>',
                 'ct': '<i>C. testosteroni</i>',
                 'ms': '<i>M. saperdae</i>',
                 'oa': '<i>O. anthropi</i>'}
        for data in fig['data']:
            try:
                data['name'] = names[data['name']]
            except KeyError:
                pass
            
        return fig

    if style == 'od_measured':
        fig.for_each_xaxis(lambda axis: axis.title.update(text='Time in hours'))
        fig.for_each_yaxis(lambda axis: axis.title.update(text='Measured OD'))

        # If file exists with sample times vlines are added
        f_times = join('data', e, 'sample_times.txt')
        if exists(f_times):
            with open(f_times, 'r') as handle:
                sample_times = handle.read().rstrip().split(',')

            for sample_time in sample_times:
                fig.add_vline(x=sample_time)

    if style == 'cfus':
        fig.for_each_xaxis(lambda axis: axis.title.update(text='Time in hours'))
        fig.for_each_yaxis(lambda axis: axis.title.update(text='CFUs/mL'))

        fig = species_colors(fig)
        fig = species_names(fig)

    if style == 'composition':
        fig.for_each_xaxis(lambda axis: axis.title.update(text='Time in hours'))
        fig.for_each_yaxis(lambda axis: axis.title.update(text='Species composition in %'))

        fig = species_colors(fig)
        fig = species_names(fig)



    fig.update_layout(font={'size': fontsize},
        title='Experiment ' + e)

    return fig


if __name__ == '__main__':
    df, fig = main()
