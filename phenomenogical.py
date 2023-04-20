from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
import plotly.express as px
import pandas as pd

rCCt = 0.3
rCOa = 0.34
rTOa = 0.32
rGOa = 0.32
rGAt = 0.27
KCOa = 7027
KGOa = 7027
KCCt = 6874
KGAt = 6500
KTOa = 0.015
qCCt = 909280 / 4
qCOa = 4321535 / 4
qGOa = 4321535 / 8
qGAt = 4321535 / 4
qTOa = 1E15
MC = 10000
MG = 15000
xs = np.arange(0, 2000, 0.5)
D = 0.15


def thiamin(y, t):
    C = y[0]
    G = y[1]
    T = y[2]
    Ct = y[3]
    Oa = y[4]
    At = y[5]
    J_COa = rCOa * C / (KCOa + C)
    J_GOa = rGOa * G / (KGOa + G)
    J_CCt = rCCt * C / (KCCt + C)
    J_GAt = rGAt * G / (KGAt + G)
    J_TOa = rTOa * T / (KTOa + T)
    dC = D * (MC - C) - J_CCt * Ct / qCCt - J_COa * Oa / qCOa
    dG = D * (MG - G) - J_GAt * At / qGAt - J_GOa * Oa / qGOa
    dT = 5E-6 * (J_GAt * At / qGAt + J_CCt * Ct / qCCt) - \
        D * T - J_TOa * Oa / qTOa
    dCt = J_CCt * Ct - D * Ct
    dAt = J_GAt * At - D * At
    dOa = min([J_TOa, max([J_GOa, J_COa])]) * Oa - D * Oa
    return [dC, dG, dT, dCt, dOa, dAt]


y = odeint(thiamin, [MC, MG, 10, 1E9, 10E9, 10E9], xs)
C, G, T, Ct, Oa, At = y[:, 0], y[:, 1], y[:, 2], y[:, 3], y[:, 4], y[:, 5]
i = len(y)
df = pd.DataFrame()
df['x'], df['y'], df['species'] = np.concatenate([xs, xs, xs, xs]), np.concatenate([Ct, Oa, At, Ct + Oa + At]),\
    np.concatenate([i * ['Ct'], i*['Oa'], i*['At'], i*['total']])

fig = px.line(df, x='x', y='y', color='species', log_y=True)
fig.show()
