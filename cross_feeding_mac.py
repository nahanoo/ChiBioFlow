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
qOa = 7.24E-5 * 0.8
qCt = 7.24E-5
KOa = 14300
KCt = 9300
NOa = 0.08
NCt = 0.08
D = 0.11
RS = 150 * 50

# Calucalted paramters
ROa = NOa / qOa
RCt = NCt / qCt 
KROa = RCt * uMax * RS / (0.1 * (KOa + RS)) - RCt
KRCt = ROa * uMax * RS / (0.1 * (KCt + RS)) - ROa


def model(y, t):
    R = y[0]
    RCt = y[1]
    Oa = y[2]
    Ct = y[3]
    dOa = uMax * R / (KOa + R) * RCt / (KROa + RCt) * Oa - D * Oa
    dCt = uCt * Ct * R / (KCt + R) - D * Ct
    dR = D * M - D * R - Oa / qOa * uOa * R / \
        (KOa + R) - Ct / qCt * uCt * R / (KCt + R)
    dRCt = - uMax * RCt / (KROa + RCt) * Oa / qOa - D * RCt + Ct / qCt
    return [dR, dRCt, dOa, dCt]


x = range(10000)
y = odeint(model, [M,NCt / qCt, NOa, NCt], x)
ODOa = y[:, 2]
ODCt = y[:, 3]


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
#plot_R()
