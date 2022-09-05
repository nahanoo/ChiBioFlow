import math
import pandas as pd
import plotly.express as px
import numpy as np
import pandas as pd
from os.path import join


def growth_rate(t0, t1, n0, n1):
    try:
        return ((math.log10(n1) - math.log10(n0)) * 2.303) / (3600*(t1-t0))
    except ValueError:
        return 0

def growth_function(t0, t1, n0, r):
    return n0 * (1 + r)**(t1 - t0)


def time_array(t0, t1, steps):
    return np.arange(t0, t1, 1/steps)


def plot(xs, ys):
    px.line(x=xs, y=ys).show()


def dilution(n1, dilution_rate, v):
    dilution_factor = (v + dilution_rate) / v
    return n1 / dilution_factor


def plot_dilutions(t0, n0, r, cycles, dilution_rate):
    x = np.ndarray(0)
    y = np.ndarray(0)

    for cycle in range(1, cycles + 1):
        t1 = cycle
        xs = time_array(t0, t1, 60)
        ys = [growth_function(t0, x, n0, r) for x in xs]
        x = np.concatenate([x, xs])
        y = np.concatenate([y, ys])
        t0 = t1
        n0 = dilution(ys[-1], dilution_rate, 20)
    plot(x,y)
    return x, y


def plot_exp_rates(csv):
    df = pd.read_csv(csv)[['exp_time', 'od_measured']]
    df["exp_time"] = df["exp_time"] / 60 / 60
    start = 0
    xs = []
    ys = []
    for i in df.index:
        if (i%60 == 0) & (i !=0):
            t0 = df.loc[start]['exp_time']
            t1 = df.loc[i-1]['exp_time']
            n0 = df.loc[start]['od_measured']
            n1 = df.loc[i-1]['od_measured']
            r = growth_rate(t0, t1, n0, n1)
            print(t0,t1,n0,n1,r)    

            xs.append(t1)
            ys.append(growth_rate(t0, t1, n0, n1))
            start = i
    plot(xs,ys)       
    return df



dilution_rate = 6.28
t0 = 0
n0 = 0.08
r = 0.32
x, y = plot_dilutions(t0,n0,r,20,dilution_rate)
plot(x,y)

