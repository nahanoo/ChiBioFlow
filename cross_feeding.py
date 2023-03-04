from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import plotly.express as px

# Measured values
M = 191 * 50
uMax = 0.32
uOa = 0.36
uCt = 0.24
qOa = 7.24E-5
qCt = 7.24E-5
KOa = 8000
KCt = 6000
NOa = 0.08
NCt = 0.08
D = 0.13
RS = 180 * 50

# Calucalted paramters
ROa = NOa / qOa
RCt = NCt / qCt 
KROa = RCt * uMax * RS / (0.15 * (KOa + RS)) - RCt
KRCt = ROa * uMax * RS / (0.15 * (KCt + RS)) - ROa


def model(y, t):
    R = y[0]
    Oa = y[1]
    Ct = y[2]
    ROa = Oa / qOa
    RCt = Ct / qCt
    dOa = uMax * R / (KOa + R) * RCt / (KROa + RCt) * Oa - D * Oa
    dCt = uMax * R / (KCt + R) * ROa / (KRCt + ROa) * Ct - D * Ct
    dR = D * M - D * R - Oa / qOa * uOa * R / \
        (KOa + R) - Ct / qCt * uCt * R / (KCt + R)
    return [dR, dOa, dCt]


x = range(200)
y = odeint(model, [M, NOa, NCt], x)
ODOa = y[:, 1]
ODCt = y[:, 2]


def plot_Com():
    plt.plot(x, ODOa, label='Oa')
    plt.plot(x, ODCt, label='Ct')
    plt.plot(x, ODCt + ODOa, label='Com')
    plt.legend()
    plt.show()


def plot_R():
    plt.plot(x, y[:, 0], label='Citrate')
    plt.legend()
    plt.show()


plot_Com()
# plot_R()
