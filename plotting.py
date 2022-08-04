from distutils.log import Log
from genericpath import exists
import plotly.express as px
import plotly.graph_objects as go
import argparse
from os.path import join
import pandas as pd
from glob import glob
from statistics import stdev, mean
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


def style_plot(fig, style, fontsize=20):

    def species_colors(fig):
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
        fig.update_layout(font={'size': fontsize},
                          xaxis_title='Time in hours',
                          yaxis_title='Measured OD')
        f_times = join('data', e, 'sample_times.txt')
        if exists(f_times):
            with open(f_times, 'r') as handle:
                sample_times = handle.read().rstrip().split(',')

            for sample_time in sample_times:
                fig.add_vline(x=sample_time)

    if style == 'cfus':
        fig.update_layout(font={'size': fontsize},
                          xaxis_title='Time in hours',
                          yaxis_title='CFUs/mL')
        fig = species_colors(fig)
        fig = species_names(fig)

    return fig


def plot_chibio():
    """Creates lineplot for parsed parameter e.g. od_measured.
    Plots every reactor as subplot. CSVs can also be parsed using
    the optional --csv flag.
    """
    df = pd.DataFrame(columns=["exp_time", "reactor", c])

    for reactor in chain:
        f = glob(join("data", e, reactor, "*data.csv"))[0]
        data = pd.read_csv(f, usecols=["exp_time", c])
        data.insert(1, "reactor", reactor)
        # time is in seconds, deviding by 60**2 to get hours
        data["exp_time"] = data["exp_time"] / 60 / 60
        df = pd.concat([df, data])
    fig = px.line(df, x="exp_time", y=c, facet_col="reactor",
                  facet_col_wrap=2, hover_data=['exp_time'])

    return fig


def plot_species(df):
    fig = px.line(df, x="sample_time", y='average', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, error_y='stdev', color='species', log_y=True)
    return fig


def plot_composition(df):
    fig = px.line(df, x="sample_time", y='composition', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, color='species',log_y=True)
    return fig


def cfu_parser():
    fs = glob(join('data', e, 'cfus*.xlsx'))
    df = pd.DataFrame(
        columns=['species', 'reactor', 'sample_time', 'dilution', 'count', 'comment'])
    for f in fs:
        species = f[-7:-5]
        cfus = pd.read_excel(f, header=1)
        cfus.insert(0, 'species', species)
        df = pd.concat([df, cfus])
    df.index = range(len(df))
    counts = []
    for count, dilution in zip(df['count'], df['dilution']):
        counts.append([int(n)/10**dilution * 0.5E2 for n in count.split('|')])
    df['average'] = [mean(count) for count in counts]
    df['stdev'] = [stdev(count) for count in counts]
    df.insert(len(df.columns), 'total', None)
    for t in df['sample_time']:
        for reactor in chain:
            tmp = df[df['sample_time'] == t]
            tmp = tmp[tmp['reactor'] == reactor]
            total = tmp['average'].sum()
            for i in tmp.index:
                df.at[i, 'total'] = total
    df.insert(len(df.columns), 'composition', None)
    df['composition'] = 100 / df['total'] * df['average']
    df = df.astype({'total': 'float64', 'composition': 'float64'})

    return df


args = parse_args()
e = args.experiment
c = args.column
mode = args.mode
with open(join('data', e, 'order.txt'), 'r') as f:
    chain = f.read().rstrip().split(',')

if mode == 'chibio':
    fig = plot_chibio()
    fig = style_plot(fig, 'od_measured')
    # fig.show()

if mode == 'species':
    df = cfu_parser()
    fig = plot_species(df)
    fig = style_plot(fig, 'cfus')
    # fig.show()

if mode == 'composition':
    df = cfu_parser()
    fig = plot_composition(df)
    fig = style_plot(fig, 'cfus')
