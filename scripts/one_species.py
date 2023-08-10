from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import plotly.express as px

"""r = 0.28082191780821913
K = 31737.470614920447
q = 90928"""
N0 = 1E7
D = 0.04
M = 10000

r = 0.42549668874172186
K = 52774
q = 432153.902910041

"""r = 0.3
q = 0.075
K = 5
M = 10
D = 0.15
N0 = 0.6
"""

def model(y, t):
    R = y[0]
    N1 = y[1]
    dR = D * M - D * R - N1 / q * r * R / \
        (R + K)
    dN1 = r * R / (K + R) * N1 - D * N1
    return [dR, dN1]


xs = np.arange(0, 4000, 1)
Y = odeint(model, [M, N0], xs)


def plot_R():
    plt.plot(xs, Y[:, 0])
    plt.show()


def plot_N():
    px.line(x=xs, y=Y[:, 1]).show()

plot_N()
