from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import plotly.express as px
from os.path import join

rCt = 0.28
rOa = 0.425
KCOa = 7000
KCCt = 7000
qCt = 90928
qOa = 432153
M = 10000
KT = 0.015
N0 = [1E8, 1E8]
xs = np.arange(0, 500, 0.5)


def plot_N(ct, oa):
    plt.plot(xs, ct, label='ct')
    plt.plot(xs, oa, label='oa')
    # plt.yscale("log")
    plt.legend()
    plt.show()


def plot_rel(ct, oa):
    plt.plot(xs, ct, label='ct')
    plt.plot(xs, oa, label='oa')
    # plt.yscale("log")
    plt.legend()
    plt.show()


def plot_R(R, T):
    plt.plot(xs, R, label='Citrate')
    plt.plot(xs, T, label='Thiamin')
    plt.yscale("log")
    plt.legend()
    plt.show()


def fixed_thiamin(y, t, c, a, D):
    R = y[0]
    T = y[1]
    Ct = y[2]
    Oa = y[3]
    dR = D * M - D * R - rCt * R / (R + KCCt) * \
        Ct/qCt - rOa * R / (R + KCOa) * Oa/qOa
    dCt = rCt * R / (KCCt + R) * Ct - D * Ct
    dOa = min(rOa * R / (KCOa + R), rOa * T / (T + KT)) * Oa - D * Oa
    #dOa = rOa * R / (KCOa + R)* rOa * T / (T + KT) * Oa - D * Oa
    if c is not None:
        dT = D * (c - T)
    if a is not None:
        dT = rCt * R / (R + KCCt) * Ct/qCt * a - D * \
            T - rOa * T / (T + KT) * Oa / 1E15
    return [dR, dT, dCt, dOa]

def ct(y,t,D):
    R = y[0]
    Ct = y[1]
    dCt = rCt * R / (KCCt + R) * Ct - D * Ct
    dR = D * M - D * R - rCt * R / (R + KCCt) * Ct/qCt
    return [dR,dCt]

def oa(y,t,D,c):
    R = y[0]
    Oa = y[1]
    T = c
    dOa = min(rOa * R / (KCOa + R), rOa * T / (T + KT)) * Oa - D * Oa
    dR = D * M - D * R - rOa * R / (R + KCOa) * Oa/qOa
    dT = D * c - D * T - rOa * T / (T + KT) * Oa / 1E15
    return [dR,dOa]


def constant_thiamine(D):
    colors = {'<i>O. anthropi</i>': '#e27e50',
              '<i>C. testosteroni</i>': '#8872cd'}
    cs = [1.5]
    dfs = []
    xs = range(800)
    for i, c in enumerate(cs):
        y0 = [M, 0, 1E8,1E8]
        y = odeint(fixed_thiamin, y0, xs, args=(c, None, D))
        df = pd.DataFrame(columns=['x', 'y', 'species', 'linegroup'])
        df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,
                                                                                   2], '<i>C. testosteroni</i>', c, y[:, 1], y[:, 0]
        dfs.append(df)
        df = pd.DataFrame(columns=['x', 'y', 'species', 'linegroup'])
        df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,
                                                                                   3], '<i>O. anthropi</i>', c, y[:, 1], y[:, 0]
        dfs.append(df)
        df = pd.DataFrame(columns=['x', 'y', 'species', 'linegroup'])
        df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,
                                                                               3] + y[:, 2], 'total', None, y[:, 1], y[:, 0]
        dfs.append(df)
    out = pd.concat(dfs)
    fig = px.line(out, x='x', y='y',
                  color='species', log_y=True,height=400,width=500)
    for i, d in enumerate(fig['data'][:-1]):
        d['line']['color'] = colors[d['legendgroup']]
    fig.update_xaxes(title='Time [h]')
    fig.update_yaxes(title='CFUs/mL')
    
    # fig.show()
    return fig

def thiamine_cross_feeding(D,xs):
    colors = {'<i>O. anthropi</i> com': '#e27e50',
              '<i>O. anthropi</i> mono': '#f25005',
              '<i>C. testosteroni</i> com': '#8872cd',
              '<i>C. testosteroni</i> mono': '#4725b0'}
    dfs = []
    y0 = [M, 0, 1E8,1E8]
    y = odeint(fixed_thiamin, y0, xs, args=(None, 5E-6, D))
    df = pd.DataFrame()
    df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,
                                                                               2], '<i>C. testosteroni</i> com', None, y[:, 1], y[:, 0]
    dfs.append(df)

    df = pd.DataFrame()
    df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,
                                                                               3], '<i>O. anthropi</i> com', None, y[:, 1], y[:, 0]
    dfs.append(df)

    """df = pd.DataFrame(columns=['x', 'y', 'species', 'linegroup'])
    df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,
                                                                               3] + y[:, 2], 'total', None, y[:, 1], y[:, 0]"""
    y0 = [M,1E8]
    y = odeint(ct,y0,xs,args=(D,))
    df = pd.DataFrame()
    df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,1], '<i>C. testosteroni</i> mono', None, None, y[:, 0]
    dfs.append(df)

    y0 = [M,1E8]
    y = odeint(oa,y0,xs,args=(D,1.5))
    df = pd.DataFrame()
    df['x'], df['y'], df['species'], df['linegroup'], df['T'], df['C'] = xs, y[:,1], '<i>O. anthropi</i> mono', None, None, y[:, 0]
    dfs.append(df)



    out = pd.concat(dfs)
    fig = px.line(out, x='x', y='y',
                  color='species', log_y=True,height=400,width=500)
    for i, d in enumerate(fig['data']):
        d['line']['color'] = colors[d['legendgroup']]
    fig.update_xaxes(title='Time [h]')
    fig.update_yaxes(title='CFUs/mL')

    fig_T = px.line(out[['x', 'T', 'linegroup']
                        ].drop_duplicates(), x='x', y='T', log_y=True)
    fig_C = px.line(out[['x', 'C', 'linegroup']
                        ].drop_duplicates(), x='x', y='C')
    return fig

