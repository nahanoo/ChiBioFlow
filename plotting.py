from genericpath import exists
import plotly.express as px
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
        fig = plot_chibio()
        fig = style_plot(fig, 'od_measured')
        fig.show()

    if mode == 'species':
        df = cfu_parser()
        fig = plot_species(df)
        fig = style_plot(fig, 'cfus')
        fig.show()

    if mode == 'composition':
        df = cfu_parser()
        fig = plot_composition(df)
        fig = style_plot(fig, 'cfus')
        fig.show()

    if mode == 'total':
        df = cfu_parser()
        fig = plot_total(df)
        #fig = style_plot(fig, 'cfus')
        fig.show()
    return df, fig


def plot_chibio():
    """Creates lineplot for parsed column e.g. od_measured.
    Plots every reactor as subplot.
    """
    df = pd.DataFrame(columns=["exp_time", "reactor", c])

    for reactor in chain:
        f = glob(join("data", e, reactor, "*data.csv"))[0]
        data = pd.read_csv(f, usecols=["exp_time", c])
        data.insert(1, "reactor", reactor)
        data["exp_time"] = data["exp_time"] / 60 / 60
        df = pd.concat([df, data])

    fig = px.line(df, x="exp_time", y=c, facet_col="reactor",
                  facet_col_wrap=2, hover_data=['exp_time'])

    return fig

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
    df['average'] = df['average'].replace(0,np.nan)
    fig = px.line(df, x="sample_time", y='average', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, error_y='stdev', color='species', log_y=True)
    print(df)
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
        fig.update_layout(font={'size': fontsize},
                          xaxis_title='Time in hours',
                          yaxis_title='Measured OD')

        # If file exists with sample times vlines are added
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

    fig.update_layout(title='Experiment ' + e)

    return fig


def cfu_parser():
    """Parses counted CFUs based on template xlsx"""
    # For eveery species there is one xlsx sheet
    fs = glob(join('data', e, 'cfus*.xlsx'))
    # df for concatenating species sheets
    df = pd.DataFrame(
        columns=['species', 'reactor', 'sample_time', 'dilution', 'count', 'comment'])
    for f in fs:
        # Species name
        species = f[-7:-5]
        cfus = pd.read_excel(f, header=1)
        cfus.insert(0, 'species', species)
        df = pd.concat([df, cfus])
    # Reindexing because of concat
    df.index = range(len(df))

    # Splitting "|" separated triplicates
    counts = []
    for count, dilution in zip(df['count'], df['dilution']):
        # Conversion to CFUs/mL sample volume 5 uL
        counts.append([int(n)/10**dilution * 0.5E2 for n in count.split('|')])
    # Calculating average and stdev
    df['average'] = [mean(count) for count in counts]
    df['stdev'] = [stdev(count) for count in counts]

    # Adding total counts for composition
    df.insert(len(df.columns), 'total', None)
    for t in df['sample_time']:
        # Subsetting for reactor and sample time
        for reactor in chain:
            tmp = df[df['sample_time'] == t]
            tmp = tmp[tmp['reactor'] == reactor]
            total = tmp['average'].sum()
            for i in tmp.index:
                df.at[i, 'total'] = total
    # Adding composition
    df.insert(len(df.columns), 'composition', None)
    df['composition'] = 100 / df['total'] * df['average']

    return df


if __name__ == '__main__':
    df,fig = main()
