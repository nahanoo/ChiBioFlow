from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import plotly.express as px

r = np.array([0.24, 0.38])
K = np.array([10962, 2000])
q = np.array([7.387404930499518e-05, 7.387404930499518e-05])
M = 11553.3127304115
D = 0.1
labels = ['Ct','Oa']


def model(y, t):
    R = y[0]
    Ns = y[1:]
    dRdN = [None for i in range(len(y))]
    DMDR = D * M - D * R
    for i in range(len(Ns)):
        DMDR = DMDR - Ns[i] / q[i] * r[i] * R / (R + K[i])
        dRdN[i + 1] = r[i] * R / (K[i] + R) * Ns[i] - D * Ns[i]
    dRdN[0] = DMDR
    return dRdN

xs = np.arange(0,150,0.5)
Y = odeint(model,[M,0.1,0.1],xs)

def plot_R():
    plt.plot(xs,Y[:,0])
    plt.show()

def plot_N():
    plt.plot(xs,Y[:,1],label=labels[0])
    plt.plot(xs,Y[:,2],label=labels[1])
    plt.legend()
    plt.show()

#plot_R()
plot_N()