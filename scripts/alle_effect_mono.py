from scipy.integrate import odeint
import numpy as np
import pandas as pd
import plotly.express as px
from os.path import join

# Two resources species can only grow on one of them
# Vector [R1_Ct,R2_Oa]
r = [0.28,0.425]

# [R1_Ct,R2_Oa,T_Oa]
K = [7, 7, 1.5E-5]

# [R1_Ct,R2_Oa,T_Oa]
q = [90928000, 432153000, 1E18]

# [R1,R2,T]
M = [10, 0]

# Cell densities
N0 = [1E8, 1E8]

u = [0, 0]

# Production rate
a = 5E-6

# Simulation parameters
xs = np.arange(0, 2000, 0.5)
y0 = M + N0


def model(y, t, D):
    R1, T, Ct, Oa = y
    uCt = r[0] * R1 / (R1 + K[0])
    uOa = min(r[1] * R1 / (R1 + K[1]), r[1] * T / (T + K[2]))

    dR1 = D * M[0] - D * R1 - uCt * Ct / q[0] - uOa * Oa / q[1]

    dOa = Oa * uOa - D * Oa
    dCt = Ct * uCt - D * Ct

    dT = a * uCt * Ct / q[0] - D * T - Oa/q[2] * r[1] * T / (T + K[2])
    return [dR1, dT, dCt, dOa]


def rates(y,t,D):
    f = model(y,t,D)
    return(f[2],f[3])

def rates_chem(y,t,D):
    f = model(y,t,D)
    return(f[3],f[4])


def simulations():
    Ds = np.arange(0.0, 0.18, 0.02)
    for D in Ds:
        y = odeint(model, y0, xs, args=(D,))

        Oa = pd.DataFrame()
        Oa['N'], Oa['species'], Oa['x'] = y[:, 3], 'Oa', xs

        Ct = pd.DataFrame()
        Ct['N'], Ct['species'], Ct['x'] = y[:, 2], 'Ct', xs

        df = pd.concat([Ct, Oa])
        fig = px.line(df, x='x', y='N', color='species')
        fig.show()

def simulation(D):
    y = odeint(model, y0, xs, args=(D,))

    Oa = pd.DataFrame()
    Oa['N'], Oa['species'], Oa['x'] = y[:, 3], 'Oa', xs

    Ct = pd.DataFrame()
    Ct['N'], Ct['species'], Ct['x'] = y[:, 2], 'Ct', xs

    df = pd.concat([Ct, Oa])
    fig = px.line(df, x='x', y='N', color='species')
    fig.show()

def alle_effect():
    df = pd.DataFrame(columns=['N','u','species'])
    D = 0
    y = odeint(model, y0, xs, args=(D,))
    for j,i in enumerate(y):
        N = i[3] + i[2]
        uCt,uOa = rates(i,xs[j],D)
        df.loc[len(df)] = [N,uCt/i[2],'Ct']
        df.loc[len(df)] = [N,uOa/i[3],'Oa']

    fig = px.line(df,x='N',y='u',color='species')
    fig.show()

def alle_effect_chemostat(D):
    df = pd.DataFrame(columns=['N','u','species'])
    y = odeint(model, y0, xs, args=(D,))
    for j,i in enumerate(y):
        N = i[3] + i[4]
        dCt,dOa = rates_chem(i,xs[j],D)
        df.loc[len(df)] = [N,dCt/i[3],'Ct']
        df.loc[len(df)] = [N,dOa/i[4],'Oa']

    fig = px.line(df,x='N',y='u',color='species')
    fig.show()

