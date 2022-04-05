import plotly.express as px
import argparse
from os.path import join, split
import pandas as pd
from glob import glob

colors = {'at': '#2c8c5a',
          'ct': '#8872cd',
          'oa': '#e27e50',
          'ms': '#e5b008'}


def parse_args():
    """Parsing variables for plotting."""
    parser = argparse.ArgumentParser(
        description="Plotting library for ChiBio.")
    parser.add_argument("experiment", help="name of the experiment directory")
    parser.add_argument("mode", help="use 'chibio' for plotting a column of the data files; \
        use 'strain' to plot strain composition.")
    parser.add_argument(
        "--column", help="column name to plot from ChiBio csv data")
    parser.add_argument(
        "--csv",
        help="path to csv if specific csv should be plotted.",
    )
    return parser.parse_args()


def plot_chibio(csv=None, transfers=False):
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
        fig = px.line(df, x="exp_time", y=c, facet_col="reactor", facet_col_wrap=2,
                      category_orders={'reactor': sorted(reactors)})
    else:
        print(csv)
        df = pd.read_csv(csv, usecols=["exp_time", c])
        df["exp_time"] = df["exp_time"] / 60 / 60
        fig = px.line(df, x="exp_time", y=c)

    if transfers:
        target_od = 0.12
        for t, od in zip(df['exp_time'], df['od_measured']):
            if od > target_od:
                fig.add_vline(x=t)
    fig.show()


def average_cfus(series):
    conversion = 20 * 1E2
    dilution = [1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7]
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
                  category_orders={'reactor': sorted(reactors)}, log_y=log, color_discrete_map=colors)
    fig.show()
    f = join('/home', 'eric', 'notes', 'talks',
             'labmeetin_2022_04_13', 'pictures', 'strains.png')
    #f = join('/home', 'eric', 'notes', 'talks',
    #         'labmeetin_2022_04_13', 'pictures', 'biofilm_strains.png')
    fig.write_image(f, scale=2)


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
        fig.show()
        f = join('/home', 'eric', 'notes', 'talks', 'labmeetin_2022_04_13',
                 'pictures', 'growth_curve_' + strain + '.png')
        fig.write_image(f, scale=2)


args = parse_args()
e = args.experiment
c = args.column
mode = args.mode
csv = args.csv
if mode == 'chibio':
    if csv is not None:
        plot_chibio(csv)
    else:
        plot_chibio()
if mode == 'strains':
    plot_strains()
    plot_strains(log=False)
if mode == 'strain':
    plot_strain()
    # plot_strains(log=False)
