import math
import pandas as pd
import plotly.express as px
import numpy as np


def growth_rate(t0, t1, n0, n1):
    return = (math.log10(t1) - math.log10(t0) * 2.303)\
        / (t1-t0)


def generation_time(r):
    return (math.log(2))/math.log(1+r)


def growth_function(t0, t1, n0, n1, r):
    return n0 * (1 + r)**(t1 - t0)


def time_array(t0, t1, steps):
    return np.arange(t0, t1, steps)

def plot(xs,ys):
    px.line(x=xs,y=ys).show()


def main():
    t0 = 0
    t1 = 1
    steps = 60
    n0 = 0.15
    n1 = 0.2

    r = growth_rate(t0, t1, n0, n1)
    xs = time_array(t0,t1,steps)
    ys = [growth_function(x) for x in xs]

    plot(xs,ys)


main()

"""28 Celsius paramters
r = growth_rate(0.93, 0.159, 0.205)
g = generation_time(r)
"""

