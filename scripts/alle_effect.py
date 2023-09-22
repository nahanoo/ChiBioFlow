from scipy.integrate import odeint
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from math import log

rCt = 0.24
rOa = 0.43
KCt = 7
KOa = 7
qCt = 0.065
qOa = 0.11
M = 10
T0 = 0
T_max = 0.00148
KT = T_max / 2
Ct0 = 0.05
Oa0 = 0.05
a = 2 * T_max / M


def model(y, t, D):
    C, T, Ct, Oa = y[0], y[1], y[2], y[3]
    J_Oa = C / (C + KOa) * T / (KT + T)
    dC = D * M - rCt * C / (C + KCt) * Ct/qCt - rOa * J_Oa * Oa / qOa - D * C
    dCt = rCt * C / (C + KCt) * Ct - D * Ct
    dOa = rOa * J_Oa * Oa - D * Oa
    dT = a * rCt * C / (C + KCt) * Ct / qCt - rOa * \
        J_Oa * Oa / (1/T_max) - D * T
    return dC, dT, dCt, dOa

def mac_arthur(y,t,D):
    C, Ct, Oa = y[0], y[1], y[2]
    J_Oa = C / (C + KOa)
    dC = D * M - rCt * C / (C + KCt) * Ct/qCt - rOa * J_Oa * Oa / qOa - D * C
    dCt = rCt * C / (C + KCt) * Ct - D * Ct
    dOa = rOa * J_Oa * Oa - D * Oa
    return dC,dCt, dOa


def growth_rates(y, t, D):
    C, T, Ct, Oa = y[0], y[1], y[2], y[3]
    uCt = rCt * C / (C + KCt)
    J_Oa = C / (C + KOa) * T / (KT + T)
    uOa = rOa * J_Oa
    return uCt, uOa


def simulation(D,t1):
    xs = np.linspace(0, t1, t1*10)
    y = odeint(model, [M, T0, Ct0, Oa0], xs, args=(D,))
    C, T, Ct, Oa = y[:, 0], y[:, 1], y[:, 2], y[:, 3]
    dfs = []
    df = pd.DataFrame(columns=['x','y','species'])
    df['x'],df['y'],df['species'] = xs,Ct,'ct'
    dfs.append(df)
    df = pd.DataFrame(columns=['x','y','species'])
    df['x'],df['y'],df['species'] = xs,Oa,'oa'
    dfs.append(df)
    fig = px.line(pd.concat(dfs),x='x',y='y',color='species')
    return fig,y

def simulation_mac(D,t1):
    xs = np.linspace(0, t1, t1*10)
    y = odeint(mac_arthur, [M, Ct0, Oa0], xs, args=(D,))
    C, Ct, Oa = y[:, 0], y[:, 1], y[:, 2]
    dfs = []
    df = pd.DataFrame(columns=['x','y','species'])
    df['x'],df['y'],df['species'] = xs,Ct,'ct'
    dfs.append(df)
    df = pd.DataFrame(columns=['x','y','species'])
    df['x'],df['y'],df['species'] = xs,Oa,'oa'
    dfs.append(df)
    fig = px.line(pd.concat(dfs),x='x',y='y',color='species')
    return fig




def alle_effect():
    df = pd.DataFrame(columns=['N', 'u', 'species'])
    D = 0.0
    t1 = 100
    xs = np.linspace(0, t1, t1*10)
    y = odeint(model, [M, T0, 0.01, 0.01], xs, args=(D,))

    for j, i in enumerate(y):
        uCt, uOa = growth_rates(i, xs[j], D)
        df.loc[len(df)] = [i[2]+i[3], uCt, 'ct']
        df.loc[len(df)] = [i[3]+i[2], uOa, 'oa']
    fig = px.line(df, x='N', y='u', color='species')
    return fig


def wash_out(Ct, Oa, D):
    def max_rate(N0, N, t, D):
        return D + 1/t * log(N/N0)
    uCt = max_rate(Ct[-100], Ct[-1], 99, D)
    uOa = max_rate(Oa[-100], Oa[-1], 99, D)
    return uCt, uOa


def alle_threshold():
    Ds = np.linspace(0.0, 0.15, 500)
    t1 = 10000
    xs = np.linspace(0, t1, t1*10)
    df = pd.DataFrame(columns=['N', 'u', 'species', 'D'])
    for D in Ds:
        y = odeint(model, [M, T0, Ct0, Oa0], xs, args=(D,))
        C, T, Ct, Oa = y[:, 0], y[:, 1], y[:, 2], y[:, 3]
        df.loc[len(df)] = [Ct[-1] + Oa[-1], growth_rates(y[-1],xs[-1],D)[0], 'ct', D]
        df.loc[len(df)] = [Oa[-1] + Ct[-1], growth_rates(y[-1],xs[-1],D)[1], 'oa', D]
    fig = px.line(df, x='D', y='u', color='species')
    return fig

def composition():
    Ds = np.linspace(0.01, 0.14, 500)
    t1 = 2000
    xs = np.linspace(0, t1, t1*10)
    df = pd.DataFrame(columns=['D', 'comp', 'species'])
    for D in Ds:
        y = odeint(model, [M, T0, Ct0, Oa0], xs, args=(D,))
        C, T, Ct, Oa = y[:, 0], y[:, 1], y[:, 2], y[:, 3]
        if D > 0.12:
            Ct_comp = 1
            Oa_comp = 0
        else:
            Ct_comp = Ct[-1] / (Ct[-1] + Oa[-1])
            Oa_comp = Oa[-1] / (Ct[-1] + Oa[-1])
        df.loc[len(df)] = [D,Ct_comp,'ct']
        df.loc[len(df)] = [D,Oa_comp,'oa']
    fig = px.line(df, x='D', y='comp', color='species')
    return fig

