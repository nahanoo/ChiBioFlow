from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import plotly.express as px

r = np.array([0.3, 0.3])
K = np.array([6000, 6000])
q = np.array([875110.8180361838, 4490271.624641353])
M = 7705.7941241261615
D = 0.05
N0 = [1E8, 1E8]
labels = ['Ct', 'Oa']


def model(y, t):
    R = y[0]
    N1 = y[1]
    N2 = y[2]
    dR = D * M - D * R - N1 / q[0] * r[0] * R / \
        (R + K[0]) - N2 / q[1] * r[1] * R / (R + K[1])
    dN1 = r[0] * R / (K[0] + R) * N1 - D * N1
    dN2 = r[1] * R / (K[1] + R) * N2 - D * N2
    return [dR, dN1, dN2]


xs = np.arange(0, 2000, 0.5)
Y = odeint(model, [M, N0[0], N0[1]], xs)


def plot_R():
    plt.plot(xs, Y[:, 0])
    plt.show()


def plot_N():
    plt.plot(xs, Y[:, 1], label=labels[0])
    plt.plot(xs, Y[:, 2], label=labels[1])
    plt.legend()
    plt.show()

