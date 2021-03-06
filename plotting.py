import plotly.express as px
import argparse
from os.path import join, split
import pandas as pd
from glob import glob
import statistics
import numpy as np

colors = {'at': '#2c8c5a',
          'ct': '#8872cd',
          'oa': '#e27e50',
          'ms': '#e5b008'}

species = ['at','ct','ms','oa']

names = {'at': '<i>A. tumefaciens</i>',
              'ct': '<i>C. testosteroni</i>',
              'ms': '<i>M. saperdae</i>',
              'oa': '<i>O. anthropi</i>'}
sample_times = [19,43,67,90,109.5]


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


def plot_chibio(csv=None, transfers=False, sampling=True):
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
        fig.show()
        f = join('/home', 'eric', 'notes', 'talks', 'labmeetin_2022_04_13',
                 'pictures', 'growth_curve_' + strain + '.png')
        fig.write_image(f, scale=2)

def rewrite():
    conversion = 1E2
    out = pd.DataFrame(columns=['day', 'reactor', 'cfus','error','species'])
    reactors = [split(element)[-1] for element in glob(join("data", e, 'M*'))]
    
    i = 0
    days = []
    for reactor in reactors:
        for f in glob(join('data', e, reactor, 'cfu*.csv')):
            species_name = f.split('.')[0][-2:]
            df = pd.read_csv(f)
            for day in df.columns[1:]:
                days.append(int(day.split('_')[-1]))
                if day == 'day_8':
                    pass
                else:
                    for counter,entry in enumerate(df[day]):
                        if '|' in str(entry):
                            counts = [int(element) for element in entry.split('|')]
                            counts = [10**counter*conversion*cfu for cfu in counts]
                            cfus = statistics.mean(counts)
                            error = statistics.stdev(counts)
                            if cfus == 0:
                                cfus = np.nan
                                error = np.nan
                            out.at[i, 'reactor'] = reactor
                            out.at[i, 'day'] = int(day.split('_')[-1])
                            out.at[i, 'cfus'] = cfus
                            out.at[i,'error'] = error
                            out.at[i,'species'] = species_name
                            i += 1
                
    fig = px.line(out, x="day", y='cfus', facet_col="reactor", facet_col_wrap=4,
                  category_orders={'reactor': sorted(reactors)}, log_y=True,color= 'species',color_discrete_map=colors,labels={
                    'day':''
                  },error_y='error')
    
    fig.show()
    fig.write_html('cfus.html')

    composition = pd.DataFrame(columns=['day', 'reactor', 'fraction','species'])
    i = 0
    for reactor in reactors:
        r = out[out['reactor'] == reactor]
        for day in days:
            d = r[r['day'] == day]
            total = 0
            for s in species:
                for count in d[d['species'] == s]['cfus']:
                    if pd.isna(count):
                        pass
                    else:
                        total += count
            for s in species:
                for count in d[d['species'] == s]['cfus']:
                    if pd.isna(count):
                        pass
                    else:
                        composition.at[i,'reactor'] = reactor
                        composition.at[i,'day'] = day
                        composition.at[i,'fraction'] = count/total
                        composition.at[i,'species'] = s
                        i += 1
    composition = composition.drop_duplicates()
    fig = px.line(composition, x="day", y='fraction',color='species', facet_col="reactor", facet_col_wrap=4,
                  category_orders={'reactor': sorted(reactors)},color_discrete_map=colors,labels={
                    'day':''
                  })
    fig.show()
    fig.write_html('composition.html')
    return composition

def plot_strains_triplicates(log=False,composition=False):
    conversion = 1E2
    out = pd.DataFrame(columns=['day', 'reactor', 'at', 'ct', 'ms', 'oa','error'])
    reactors = [split(element)[-1] for element in glob(join("data", e, 'M*'))]
    
    i = 0
    days = []
    for reactor in reactors:
        for f in glob(join('data', e, reactor, 'cfu*.csv')):
            strain = f.split('.')[0][-2:]
            df = pd.read_csv(f)
            for day in df.columns[1:]:
                days.append(int(day.split('_')[-1]))
                if day == 'day_8':
                    pass
                else:
                    for counter,entry in enumerate(df[day]):
                        if '|' in str(entry):
                            counts = [int(element) for element in entry.split('|')]
                            counts = [10**counter*conversion*cfu for cfu in counts]
                            cfus = statistics.mean(counts)
                            error = statistics.stdev(counts)
                            if cfus == 0:
                                cfus = np.nan
                                error = np.nan
                            out.at[i, 'reactor'] = reactor
                            out.at[i, 'day'] = int(day.split('_')[-1])
                            out.at[i, strain] = cfus
                            out.at[i,'error'] = error
                            i += 1
                
    fig = px.line(out, x="day", y=['at', 'ct', 'ms', 'oa'], facet_col="reactor", facet_col_wrap=4,
                  category_orders={'reactor': sorted(reactors)}, log_y=log, color_discrete_map=colors,labels={
                    'day':''
                  },error_y='error')
    #fig.show()
    composition = pd.DataFrame(columns=['day', 'reactor', 'at', 'ct', 'ms', 'oa'])
    i = 0
    for reactor in reactors:
        r = out[out['reactor'] == reactor]
        for day in days:
            d = r[r['day'] == day]
            total = 0
            for s in species:
                for count in d[s]:
                    if pd.isna(count):
                        pass
                    else:
                        total += count
            for s in species:
                for count in d[s]:
                    if pd.isna(count):
                        pass
                    else:
                        composition.at[i,'reactor'] = reactor
                        composition.at[i,'day'] = day
                        composition.at[i,s] = count/total
                        i += 1
    composition = composition.drop_duplicates()
    fig = px.scatter(composition, x="day", y=['at', 'ct', 'ms', 'oa'], facet_col="reactor", facet_col_wrap=4,
                  category_orders={'reactor': sorted(reactors)},color_discrete_map=colors,labels={
                    'day':''
                  })
    fig.show()

    
    """
    fig.for_each_trace(lambda t: t.update(name = names[t.name],
                                      legendgroup = names[t.name],
                                      hovertemplate = t.hovertemplate.replace(t.name, names[t.name])
                                     )
                  )
    fig.update_layout(font={'size':40},
            xaxis_title='Day',
            yaxis_title='CFUs/mL')

    fig.show()
    """
    return out


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
    #plot_strains()
    plot_strains(log=True)
if mode == 'strain':
    plot_strain()
    # plot_strains(log=False)

if mode == 'biofilm':
    plot_strains(csv='biofilm_cfu*.csv')

if mode == 'strains_triplicates':
    out = plot_strains_triplicates()

if mode == 'rewrite':
    out = rewrite()