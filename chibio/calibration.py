import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go
from os.path import join, split
import csv
import pandas as pd
from glob import glob
import math
from matplotlib import pyplot as plt
from sympy import *

init_printing()

reactors = ['M' + str(j) for j in range(8)]
setup_name = 'calibrations_standard'


def exp_fit(x, m, t):
    return m * np.exp(-t * x)


def get_od(raw):
    raw_ods = [
        (12946.666666666666, 0.06),
        (10381.0, 0.11),
        (5177.0, 0.3),
        (1298.2222222222224, 0.72),
        (425.5555555555556, 1.1),
        (84.8888888888889, 2.0),
        (26.111111111, 3.7),
        (10.0, 6.0),
    ]
    if raw >= raw_ods[0][0]:
        return 0.05
    elif raw <= raw_ods[-1][0]:
        return 7
    for i, (r, od) in enumerate(raw_ods):
        if raw > r:
            xs = [r, raw_ods[i - 1][0]]
            ys = [od, raw_ods[i - 1][1]]
            m = (ys[1] - ys[0]) / (xs[1] - xs[0])
            b = ys[0] - m * xs[0]
            return m * raw + b


def dump_fit(plot=False):
    for reactor in reactors:
        fin = join('calibrations', reactor + '_raw.csv')
        xs = []
        ys = []
        df = pd.read_csv(fin, index_col=False)
        xs += df['od'].to_list() * 3
        for col in df.columns[1:]:
            ys += df[col].to_list()
        (m, t), cv = curve_fit(exp_fit, xs, ys)
        fout = join('calibrations', reactor + '_params.csv')
        df = pd.DataFrame(columns=['m', 't'], index=[0])
        df.at[0, 'm'] = m
        df.at[0, 't'] = t
        df.to_csv(fout, index=False)

        if plot:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=xs, y=ys, mode='markers'))
            # fig.add_trace(go.Scatter(x=df['od'],y=mono_exp(df['od'],m,t),mode='lines'))
            fig.add_trace(
                go.Scatter(x=np.arange(0, 5.1, 0.05),
                           y=exp_fit(np.arange(0, 5.1, 0.05), m, t),
                           mode='lines'))
            fig.update_layout(title=reactor)
            fig.show()


def exp_fit_all():
    {
        'average': {
            0: 12946.666666666666,
            1: 10381.0,
            2: 5177.0,
            3: 1298.2222222222224,
            4: 425.5555555555556,
            5: 84.8888888888889,
            6: 26.11111111111111,
            7: 10.0
        },
        'od': {
            0: 0.06,
            1: 0.11,
            2: 0.3,
            3: 0.72,
            4: 1.1,
            5: 2.0,
            6: 3.7,
            7: 6.0
        }
    }

    dfs = []
    fs = glob(join(setup_name, 'M*raw.csv'))
    for f in fs:
        name = split(f)[-1][:2]
        df = pd.read_csv(f)
        df.insert(len(df.columns), 'reactor', name)
        df.insert(len(df.columns), 'raw_average', name)
        for i, row in df.iterrows():
            raws = [row['raw_1'], row['raw_2'], row['raw_3']]
            avg = np.average(raws)
            df.at[i, 'average'] = avg
        #df = df.loc[:4]
        dfs.append(df)
    out = pd.concat(dfs)

    ods = set(out['od'])
    avg_df = pd.DataFrame(index=sorted(ods), columns=['average'])
    for o in ods:
        odf = out[out['od'] == o]
        avg_df.at[o, 'average'] = np.average(odf['average'])
    (m, t), cv = curve_fit(exp_fit, avg_df.index, avg_df['average'])
    avg_df['od'] = avg_df.index
    avg_df.index = range(len(avg_df))
    return avg_df


R, OD, m, k, b = symbols('R OD m k b')
R1, R2, OD1, OD2 = symbols('R1 R2 OD1 OD2')

# Define the equations based on two points with the new terms k and b
R, OD, k, b = symbols('R OD k b')
m = 0.05  # Given that m is determined at OD = 0.05

# Known data points for OD and R
R1, R2, OD1, OD2 = symbols('R1 R2 OD1 OD2')

# Set up the adjusted equation with m fixed
eq1 = Eq(OD1, k * log(m / R1, 10) + b)
eq2 = Eq(OD2, k * log(m / R2, 10) + b)

# Solve for k and b
solution = solve((eq1, eq2), (k, b))

a = exp_fit_all()
fig = go.Figure()
fs = sorted(
    glob('/home/eric/ChiBioFlow/chibio/calibrations_mitri/M[0-9]_raw.csv'))
for f in fs:
    df = pd.read_csv(f)
    df = df.sort_values('od')
    #df = df[df['od'] <= 0.9]
    df.index = range(len(df))
    avg = np.average(df[['raw_1', 'raw_2', 'raw_3']], axis=1)
    R0, R1 = [avg[0], avg[7]]
    ODs = df['od'].to_list()
    OD0, OD1 = [ODs[0], ODs[7]]
    b_num = float(solution[b].subs({
        'OD1': OD0,
        'OD2': OD1,
        'R1': R0,
        'R2': R1
    }))
    k_num = float(solution[k].subs({
        'OD1': OD0,
        'OD2': OD1,
        'R1': R0,
        'R2': R1
    }))
    y = k_num * np.log10(m / avg) + b_num
    fig.add_trace(go.Scatter(x=df['od'], y=y))
fig.show()
