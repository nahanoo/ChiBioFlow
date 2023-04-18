from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import plotly.express as px

r = 0.3
KS = 6000
qCt = 909280
qOa = 4321535
M = 10000
KT = 0.015
N0 = [1E8, 1E9]
xs = np.arange(0, 200, 0.5)


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


def fixed_thiamin(y, t, c):
    D = 0.15
    R = y[0]
    T = y[1]
    Ct = y[2]
    Oa = y[3]
    dR = D * M - D * R - r * R / (R + KS) * Ct/qCt - r * R / (R + KS) * Oa/qOa
    dT = r * R / (R + KS) * Ct/qCt * 2E-5 - D * T
    dCt = r * R / (KS + R) * Ct - D * Ct
    dOa = min(r * R / (KS + R), r * T / (T + KT)) * Oa - D * Oa
    return [dR, dT, dCt, dOa]


y = odeint(fixed_thiamin, [M, 0, N0[0], N0[1]], xs, args=(None,))
Ct, Oa = y[:, 2], y[:, 3]
R, T = y[:, 0], y[:, 1]
dfs = []
df = pd.DataFrame(
    columns=['time', 'y', 'species', 'R', 'Thiamine [\u03BCM]', 'total'])
df['time'], df['y'], df['species'], df['R'], df['Thiamine [\u03BCM]'], df['total'] = xs, Ct, '<i>C. testosteroni</i>', R, T, Ct+Oa
dfs.append(df)
df = pd.DataFrame(
    columns=['time', 'y', 'species', 'R', 'Thiamine [\u03BCM]', 'total'])
df['time'], df['y'], df['species'], df['R'], df['Thiamine [\u03BCM]'], df['total'] = xs, Oa, '<i>O. anthropi</i>', R, T, Ct+Oa
dfs.append(df)
out = pd.concat(dfs)
fig = px.line(out, x='time', y='y', color='species')
fig.show()
total = out['y'] / out['total']
fig = px.line(out, x='time', y=total, color='species')
fig.show()
dT = out[['Thiamine [\u03BCM]', 'time']].drop_duplicates()
fig = px.line(dT, x='time', y='Thiamine [\u03BCM]')
fig.show()


"""cs = np.linspace(0, 0.02, 10)
dfs = []
for c in cs:
    y = odeint(fixed_thiamin, [M, c, N0[0], N0[1]], xs, args=(c,))
    Ct, Oa = y[:, 2], y[:, 3]
    R, T = y[:, 0], y[:, 1]
    df = pd.DataFrame(
        columns=['time', 'y', 'species', 'R', 'Thiamine [\u03BCM]'])
    df['time'], df['y'], df['species'], df['R'], df['Thiamine [\u03BCM]'] = xs, Ct, '<i>C. testosteroni</i>', R, c
    dfs.append(df)
    df = pd.DataFrame(
        columns=['time', 'y', 'species', 'R', 'Thiamine [\u03BCM]'])
    df['time'], df['y'], df['species'], df['R'], df['Thiamine [\u03BCM]'] = xs, Oa, '<i>O. anthropi</i>', R, c
    dfs.append(df)
out = pd.concat(dfs)
fig = px.line(out, x='time', y='y', color='species',
              line_dash='Thiamine [\u03BCM]')
fig.update_xaxes(title='Time [h]')
fig.update_yaxes(title='CFUs/mL')
fig.show()"""

"""c = 188
comp = np.linspace(1E7, 1E10, 3)
dfs = []
for i, j in enumerate(comp):
    y0 = [M, c, comp[i], comp[-(i + 1)]]
    print(y0)
    y = odeint(fixed_thiamin, y0, xs, args=(c,))
    df = pd.DataFrame(columns=['x', 'y', 'species', 'linegroup'])
    df.loc[0] = [xs[-1], y[:,
                           2][-1], '<i>C. testosteroni</i>', i]
    df.loc[1] = [xs[-1], y[:,
                           3][-1], '<i>O. anthropi</i>', i]
    dfs.append(df)
out = pd.concat(dfs)
fig = px.scatter(out, x='linegroup', y='y', color='species')
fig.update_xaxes(type='category')
fig.show()"""

"""def model(y, t):
    R = y[0]
    T = y[1]
    Ct = y[2]
    Oa = y[3]
    dR = D * M - D * R - r * R / (R + KS) * (Ct/qCt + Oa/qOa)
    dCt = r * R / (KS + R) * Ct - D * Ct
    dOa = r * R / (KS + R) * T / (T + KT) * Oa - D * Oa
    dT = Ct * BV - D * T - r * T / (T + KT) * Oa / qV
    return [dR, dT, dCt, dOa]
"""
