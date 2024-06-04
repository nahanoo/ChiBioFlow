from scipy.integrate import odeint
import numpy as np
import pandas as pd
import plotly.express as px
from sympy import symbols, solve, Eq
from sympy import init_printing

init_printing()


rAt = 0.22
rOa = 0.15
KAt = 15.4
YTOa = 22808
YCOa = 0.053
YAt = 0.075
KCOa = 0.55
T_limOa = 1.4477458712604146e-05
KTOa = T_limOa / 2
M = 30
YTOa = 32377
a = 2* 8.82460653285622e-8


JOa, JAt, C, T, dC, dT, dAt, dOa, dC, dT, At, Oa, D, JOaC = symbols(
    "JOa,JAt,C,T,dC,dT,dAt,dOa,dC,dT,At,Oa,D,JOaC"
)
JOa = Eq(JOa, rOa * C / (C + KCOa) * T / (T + KTOa))
JOaC = Eq(JOaC, rOa * C / (C + KCOa))
JAt = Eq(JAt, rAt * C / (C + KAt))
dAt = Eq(dAt, JAt.rhs * At - D * At)
dAt = Eq(dOa, JOa.rhs * At - D * Oa)
dC = Eq(dC, M * D - D * C - JOa.rhs * Oa / YCOa - JAt.rhs * At / YAt)
dT = Eq(dT, a * JAt.rhs * At / YAt - JOa.rhs * Oa / YTOa)


def competition(y, t, D):
    C, At, Oa = y[0], y[1], y[2]
    JAt = rAt * C / (C + KAt)
    JOa = rOa * C / (C + KCOa)
    dC = D * M - JAt * At / YAt - JOa * Oa / YCOa - D * C
    dAt = JAt * At - D * At
    dOa = JOa * Oa - D * Oa
    return dC, dAt, dOa


def thiamine_supply(y, t, D, MT):
    C, T, At, Oa = y[0], y[1], y[2], y[3]
    JAt = rAt * C / (C + KAt)
    JOa = rOa * C / (C + KCOa) * T / (T + KTOa)
    dC = D * M - JAt * At / YAt - JOa * Oa / YCOa - D * C
    dAt = JAt * At - D * At
    dOa = JOa * Oa - D * Oa
    dT = D * MT - JOa * Oa / YTOa - D * T
    return dC, dT, dAt, dOa


def siumulate_thiamine_supply(T,D):
    t1 = 500
    xs = np.linspace(0, t1, t1 * 10)
    y = odeint(thiamine_supply, [M, 0, 0.05, 0.05], xs, args=(D, T))
    C, T, At, Oa = y[:, 0], y[:, 1], y[:, 2], y[:, 3]
    dfs = []
    df = pd.DataFrame(columns=["time", "N", "species"])
    df["time"], df["N"], df["species"] = xs, At, "At"
    dfs.append(df)
    df = pd.DataFrame(columns=["time", "N", "species"])
    df["time"], df["N"], df["species"] = xs, Oa, "Oa"
    dfs.append(df)
    fig = px.line(
        pd.concat(dfs), x="time", y="N", color="species", width=500, height=400
    )
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    return pd.concat(dfs)


def cross_feeding(y, t, D):
    C, T, At, Oa = y[0], y[1], y[2], y[3]
    JAt = rAt * C / (C + KAt)
    JOa = rOa * C / (C + KCOa) * T / (T + KTOa)
    dC = D * M - JAt * At / YAt - JOa * Oa / YCOa - D * C
    dAt = JAt * At - D * At
    dOa = JOa * Oa - D * Oa
    dT = a * JAt * At / YAt - JOa * Oa / YTOa
    return dC, dT, dAt, dOa


def simulate_competition(D, t1):
    xs = np.linspace(0, t1, t1 * 10)
    y = odeint(competition, [M, 0.05, 0.05], xs, args=(D,))
    C, At, Oa = y[:, 0], y[:, 1], y[:, 2]
    dfs = []
    df = pd.DataFrame(columns=["time", "N", "species"])
    df["time"], df["N"], df["species"] = xs, At, "At"
    dfs.append(df)
    df = pd.DataFrame(columns=["time", "N", "species"])
    df["time"], df["N"], df["species"] = xs, Oa, "Oa"
    dfs.append(df)
    fig = px.line(pd.concat(dfs), x="time", y="N", color="species")
    fig.show()
    return pd.concat(dfs)


def simulate_cross_feeding():
    t1 = 500
    D = 0.1
    xs = np.linspace(0, t1, t1 * 10)
    y = odeint(cross_feeding, [M, 2E-5, 0.05, 0.05], xs, args=(D,))
    C, T, At, Oa = y[:, 0], y[:, 1], y[:, 2], y[:, 3]
    dfs = []
    df = pd.DataFrame(columns=["time", "N", "species"])
    df["time"], df["N"], df["species"] = xs, At, "At"
    dfs.append(df)
    df = pd.DataFrame(columns=["time", "N", "species"])
    df["time"], df["N"], df["species"] = xs, Oa, "Oa"
    dfs.append(df)
    fig = px.line(
        pd.concat(dfs), x="time", y="N", color="species", width=500, height=400
    )
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    # fig.write_image('at_oa_alanine.svg')
    return pd.concat(dfs)


