from statistics import stdev, mean
from os.path import join
import pandas as pd
from glob import glob
import json
import numpy as np


def cfu_parser(e):
    with open(join('/','home','eric','ChiBioFlow', 'data', e, 'order.txt'), 'r') as f:
        order = f.read().rstrip().split(',')

    """Parses counted CFUs based on template xlsx"""
    # For eveery species there is one xlsx sheet
    fs = glob(join('/','home','eric','ChiBioFlow', 'data', e, 'cfus*.xlsx'))
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
        counts.append([(int(n)/10**dilution) * 100 for n in count.split('|')])
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
    df['composition'] = df['average'] / df['total']
    return df, order


def chibio_parser(e, down_sample=False, sample_size=10):
    def average(df):
        out = pd.DataFrame(
            columns=["exp_time", 'od_measured', 'FP1_base', 'FP1_emit1', 'FP1_emit2'])
        i = 0
        j = 0
        while True:
            try:
                sliece = df.loc[i:i+sample_size]
                out.loc[j] = [sliece.iloc[-1]['exp_time'], np.average(sliece['od_measured']), np.average(
                    sliece['FP1_base']), np.average(sliece['FP1_emit1']), np.average(sliece['FP1_emit2'])]
                i += sample_size
                j += 1
            except IndexError:
                break
        return out

    with open(join('/','home','eric','ChiBioFlow','data', e, 'order.txt'), 'r') as f:
        order = f.read().rstrip().split(',')
    dfs = []
    for reactor in order:
        f = glob(join('/','home','eric','ChiBioFlow','data', e, reactor, "*data.csv"))[0]
        data = pd.read_csv(
            f, usecols=["exp_time", 'od_measured', 'FP1_base', 'FP1_emit1', 'FP1_emit2'])
        if down_sample:
            data = average(data)
        data.insert(1, "reactor", reactor)
        data["exp_time"] = data["exp_time"] / 60 / 60
        dfs.append(data)

    df = pd.concat(dfs)
    dfs = []
    for c in ['od_measured', 'FP1_emit1', 'FP1_emit2','FP1_base']:
        tmp = pd.DataFrame(columns=['exp_time','reactor','measurement','sensor'])
        tmp['reactor'] = df['reactor']
        tmp['exp_time'] = df['exp_time']
        tmp['measurement'] = df[c]
        tmp['sensor'] = c
        dfs.append(tmp)
    df = pd.concat(dfs)



    return df, order
