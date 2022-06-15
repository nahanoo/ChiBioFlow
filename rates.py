import math
import pandas as pd
import plotly.express as px
import numpy as np


def growth_rate(t0, t1, n0, n1):
    try:
        return ((math.log10(n1) - math.log10(n0)) * 2.303) / (t1-t0)
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

def main(t0,n0,cycles,dilution_rate):
    x = np.ndarray(0)
    y = np.ndarray(0)

    for cycle in range(1, cycles + 1):
        t1 = cycle
        xs = time_array(t0, t1, 60)
        ys = [growth_function(t0, x, n0, r) for x in xs]
        x = np.concatenate([x,xs])
        y = np.concatenate([y,ys])
        t0 = t1
        n0 = dilution(ys[-1],dilution_rate,15)

    return x,y

v = 15
dilution_rate = 3.41
t0 = 59.046
t1 = 61.340
steps = 60
n0 = 0.118
n1 = 0.201
r = growth_rate(t0, t1, n0, n1)
t0 = 0
n0 = 0.06
x, y = main(t0,n0,10,dilution_rate)
plot(x,y)

"""28 Celsius paramters
r = growth_rate(0.93, 0.159, 0.205)
g = generation_time(r)
def generation_time(r):
    return (math.log(2))/math.log(1+r)
"""
