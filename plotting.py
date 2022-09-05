
from audioop import add
from genericpath import exists
import plotly.express as px
import plotly.graph_objects as go
import argparse
from os.path import join
from chibio_parser import cfu_parser
from chibio_parser import chibio_parser
import numpy as np
from model_class import Chain


def plot_od(e, multi=True, model=False, chain=None):
    """Creates lineplot for parsed column e.g. od_measured.
    Plots every reactor as subplot.
    """
    df, order = chibio_parser(e)

    if multi:
        fig = px.line(df, x="exp_time", y='od_measured', facet_col="reactor", color='temp', hover_data=[
                      'exp_time'], category_orders={'reactor': order})

    if not multi:
        fig = px.line(df, x="exp_time", y=c, color='temp',
                      hover_data=['exp_time'],)

    fig = style_plot(e, fig, 'od_measured')

    if model:
        temps = df[['reactor', 'temp']].drop_duplicates()['temp'].to_list()
        chain.experiment(int(max(df['exp_time'])))
        for counter, f in enumerate(fig['data']):
            xs = chain.chain[counter].xs
            ys = chain.chain[counter].ys
            name = str(temps[counter]) + ' model'
            fig.add_trace(go.Scatter(x=xs, y=ys, name=name, line=dict(dash='dash')),
                          row=1, col=counter + 1)
        fig = style_plot(e, fig, 'od_measured')

    return df, fig


def plot_total(e):
    """Plot sum of CFUs of all species"""
    df, order = cfu_parser(e)
    df = df[['reactor', 'sample_time', 'total']].drop_duplicates()
    df['total'] = df['total'].replace(0, np.nan)
    fig = px.line(df, x="sample_time", y='total', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': order}, log_y=True)
    fig.update_layout(font={'size': 20},
                      xaxis_title='Time in hours',
                      yaxis_title='CFUs/mL')
    return df, fig


def plot_species(e):
    """Plots CFUs based on parsed xlsx sheet"""
    df, order = cfu_parser(e)
    df['average'] = df['average'].replace(0, np.nan)
    fig = px.line(df, x="sample_time", y='average', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': order}, error_y='stdev', color='species', log_y=True)
    fig = style_plot(e, fig, 'cfus')
    return df, fig


def plot_composition(e):
    """Plots community composition in percent"""
    df, order = cfu_parser(e)
    fig = px.line(df, x="sample_time", y='composition', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': order}, color='species', log_y=False)
    fig = style_plot(e, order, fig, 'composition')

    return df, fig


def style_plot(e, fig, style, fontsize=14):
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
        temp_colors = {28.0: '#420f99',
                       33.0: '#a234eb',
                       38.0: '#eb34c0',
                       43.0: '#eb3434'}
        for data in fig['data']:
            temp = data['name'].split(' ')[0]
            data.line.color = temp_colors[float(temp)]
        fig.for_each_xaxis(
            lambda axis: axis.title.update(text='Time in hours'))
        fig.for_each_yaxis(lambda axis: axis.title.update(text='Measured OD'))

        # If file exists with sample times vlines are added
        f_times = join('data', e, 'sample_times.txt')
        if exists(f_times):
            with open(f_times, 'r') as handle:
                sample_times = handle.read().rstrip().split(',')

            for sample_time in sample_times:
                fig.add_vline(x=sample_time)
        fig.update_layout(
            height=350)

    if style == 'cfus':
        fig.for_each_xaxis(
            lambda axis: axis.title.update(text='Time in hours'))
        fig.for_each_yaxis(lambda axis: axis.title.update(text='CFUs/mL'))

        fig = species_colors(fig)
        fig = species_names(fig)

    if style == 'composition':
        fig.for_each_xaxis(
            lambda axis: axis.title.update(text='Time in hours'))
        fig.for_each_yaxis(lambda axis: axis.title.update(
            text='Species composition in %'))

        fig = species_colors(fig)
        fig = species_names(fig)

    fig.update_layout(font={'size': fontsize},
                      title=e)

    return fig


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
    args = parse_args()
    e = args.experiment
    c = args.column
    mode = args.mode
    with open(join('data', e, 'order.txt'), 'r') as f:
        order = f.read().rstrip().split(',')

    if mode == 'chibio_multi':
        df, fig = plot_od(e)
        fig.show()

    if mode == 'chibio_single':
        df, fig = plot_od(e, multi=False)
        fig.show()

    if mode == 'chibio_model':
        df, order = chibio_parser(e)
        temps = df[['reactor', 'temp']].drop_duplicates()['temp'].to_list()
        chain = Chain(temps)
        df, fig = plot_od(e, model=True, chain=chain)
        fig.show()
        return df, fig

    if mode == 'species':
        df, fig = plot_species(e)
        fig.show()
        return df, fig

    if mode == 'composition':
        df, fig = plot_composition(e)
        fig = style_plot(e, order, fig, 'composition')
        fig.show()
        return df, fig

    if mode == 'total':
        df, fig = plot_total(e)
        fig.show()
        return df, fig

    return df, fig


if __name__ == '__main__':
    df, fig = main()
