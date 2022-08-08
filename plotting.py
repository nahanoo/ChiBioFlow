from genericpath import exists
import plotly.express as px
import argparse
from os.path import join
import pandas as pd
from glob import glob
from statistics import stdev, mean


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


<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 714573da5466ad345b2ce34bffd315b045f08b68
def plot_chibio(csv=None, transfers=False, sampling=False):
    """Creates lineplot for parsed parameter e.g. od_measured.
    Plots every reactor as subplot. CSVs can also be parsed using
    the optional --csv flag.
    """
    reactors = [split(element)[-1] for element in glob(join("data", e, 'M*'))]
    df = pd.DataFrame(columns=["exp_time", "reactor", c])
    if csv is None:
        for reactor in reactors:
            f = glob(join("data", e, reactor, "2*.csv"))
            if len(f) > 1:
                print(
                    "There are multiple csv as sources. Clean direcotry first. Or use the optional --csv flag."
                )
                return
            data = pd.read_csv(f[0], usecols=["exp_time", c])
            data.insert(1, "reactor", reactor)
            # time is in seconds, deviding by 60**2 to get hours
            data["exp_time"] = data["exp_time"] / 60 / 60
            df = df.append(data)
        if c == 'od_measured':
            fig = px.line(df, x="exp_time", y=c, facet_col="reactor", facet_col_wrap=2,hover_data=['exp_time']  ,
                      category_orders={'reactor': sorted(reactors)})
            fig.update_layout(font={'size':20})
            fig.update_layout(
            xaxis_title='Time in hours',
            yaxis_title='Measured OD')
        else:
            fig = px.line(df, x="exp_time", y=c, facet_col="reactor", facet_col_wrap=2,
                        category_orders={'reactor': sorted(reactors)})
    else:
        print(csv)
        df = pd.read_csv(csv, usecols=["exp_time", c])
        df["exp_time"] = df["exp_time"] / 60 / 60
        fig = px.line(df, x="exp_time", y=c)
    
    if transfers:
        target_od = 0.2
        for t, od in zip(df['exp_time'], df['od_measured']):
            if od > target_od:
                fig.add_vline(x=t)

    if sampling:
        for sample_time,day in zip(sample_times,[1,2,3,4,5]):
            fig.add_vline(x=sample_time,annotation_text=day,line_color="orange")

    fig.write_html('tg_v1.html')
    fig.show()


def average_cfus(series):
    conversion = 1E2
    dilution = [1E0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6]
    cfus = []
    for counter, cfu in enumerate(series):
        if cfu == 'overgrown':
            pass
        else:
            cfu = int(cfu)
            if cfu >= 8:
                cfus.append(cfu * dilution[counter] * conversion)
    if len(cfus) == 0:
        avg = None
    else:
        avg = sum(cfus)/len(cfus)
    return avg


def plot_strains(log=True):
    out = pd.DataFrame(columns=['day', 'reactor', 'at', 'ct', 'ms', 'oa'])
    reactors = [split(element)[-1] for element in glob(join("data", e, 'M*'))]
    i = 0
    for reactor in reactors:
        for f in glob(join('data', e, reactor, 'cfu*.csv')):
            strain = f.split('.')[0][-2:]
            df = pd.read_csv(f)
            for day in df.columns[1:]:
                out.at[i, 'reactor'] = reactor
                out.at[i, 'day'] = day.split('_')[-1]
                out.at[i, strain] = average_cfus(df[day])
                i += 1
    fig = px.line(out, x="day", y=['at', 'ct', 'ms', 'oa'], facet_col="reactor", facet_col_wrap=4,
                  category_orders={'reactor': sorted(reactors)}, log_y=log, color_discrete_map=colors,labels={
                    'day':''
                  })
    fig.for_each_trace(lambda t: t.update(name = names[t.name],
                                      legendgroup = names[t.name],
                                      hovertemplate = t.hovertemplate.replace(t.name, names[t.name])
                                     )
                  )
    fig.update_layout(font={'size':20},
            xaxis_title='Day',
            yaxis_title='CFUs/mL')

    fig.show()
    #f = join('/home', 'eric', 'notes', 'talks',
    #         'labmeetin_2022_04_13', 'pictures', 'strains.png')
    # f = join('/home', 'eric', 'notes', 'talks',
    #         'labmeetin_2022_04_13', 'pictures', 'biofilm_strains.png')
    #fig.write_image(f, scale=2)


def plot_strain(log=True):
    out = pd.DataFrame(columns=['day', 'reactor', 'at', 'ct', 'ms', 'oa'])
    reactors = [split(element)[-1] for element in glob(join("data", e, 'M*'))]
    i = 0
    for reactor in reactors:
        for f in glob(join('data', e, reactor, 'cfu*.csv')):
            strain = f.split('.')[0][-2:]
            df = pd.read_csv(f)
            for day in df.columns[1:]:
                out.at[i, 'reactor'] = reactor
                out.at[i, 'day'] = day.split('_')[-1]
                out.at[i, strain] = average_cfus(df[day])
                i += 1
    for strain in ['at', 'ct', 'ms', 'oa']:
        fig = px.line(out, x="day", y=[strain], facet_col="reactor", facet_col_wrap=4,
                      category_orders={'reactor': sorted(reactors)}, color_discrete_map=colors, log_y=log)
=======
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
>>>>>>> a61bdd4a5dfc00c68a77de50c19f2ce6b8e2425e
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


def plot_species(df):
    """Plots CFUs based on parsed xlsx sheet"""
    fig = px.line(df, x="sample_time", y='average', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, error_y='stdev', color='species', log_y=True)
    return fig


def plot_composition(df):
    """Plots community composition in percent"""
    fig = px.line(df, x="sample_time", y='composition', facet_col="reactor",
                  facet_col_wrap=2, category_orders={'reactor': chain}, color='species', log_y=True)
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
    main()
