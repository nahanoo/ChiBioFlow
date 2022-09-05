from statistics import stdev, mean
from os.path import join
import pandas as pd
from glob import glob
import json

def cfu_parser(e):
    with open(join('data', e, 'order.txt'), 'r') as f:
        order = f.read().rstrip().split(',')

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
        for reactor in order:
            tmp = df[df['sample_time'] == t]
            tmp = tmp[tmp['reactor'] == reactor]
            total = tmp['average'].sum()
            for i in tmp.index:
                df.at[i, 'total'] = total
    # Adding composition
    df.insert(len(df.columns), 'composition', None)
    df['composition'] = 100 / df['total'] * df['average']

    return df, order

def chibio_parser(e,c='od_measured'):
    df = pd.DataFrame(columns=["exp_time", "reactor", c])

    with open(join('data', e, 'order.txt'), 'r') as f:
        order = f.read().rstrip().split(',')

    for reactor in order:
        f = glob(join("data", e, reactor, "*data.csv"))[0]
        data = pd.read_csv(f, usecols=["exp_time", c])
        data.insert(1, "reactor", reactor)
        data["exp_time"] = data["exp_time"] / 60 / 60
        data.insert(len(data.columns),'temp',None)
        f = glob(join("data", e, reactor, "*.txt"))[0]
        j = open(f)
        j_data = json.load(j)
        data['temp'] = j_data['Thermostat']['last']

        df = pd.concat([df, data])
    if e == 'overnight_gradient_06_16_failed':
        df = df[df['od_measured'] < 0.8]

    return df,order