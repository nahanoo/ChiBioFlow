from scipy.integrate import odeint
import numpy as np
import pandas as pd
import plotly.express as px
from os.path import join

# Two resources species can only grow on one of them
# Vector [R1_Ct,R2_Oa]
r = [0.28082191780821913,0.42549668874172186]

# [R1_Ct,R2_Oa,T_Oa]
K = [7, 7, 1.5E-5]

# [R1_Ct,R2_Oa,T_Oa]
q = [909280003.9898741, 4321534.902910041, 1E15]

# [R1,R2,T]
M = [15, 10, 0]

# Cell densities
N0 = [1E7, 1E7]

u = [0, 0]

# Production rate
a = 1E-6

# Simulation parameters
xs = np.linspace(0,200,200)
y0 = M + N0 + u


def model(y, t, D):
    R1, R2, T, Ct, Oa, uCt, uOa = y
    uCt = r[0] * R1 / (R1 + K[0])
    uOa = min((r[1] * R2 / (R2 / K[1])), r[1] * T / (T + K[2]))

    dR1 = D * M[0] - D * R1 - uCt * Ct / q[0]
    dR2 = D * M[1] - D * R2 - uOa * Oa / q[1]

    dOa = Oa * uOa - D * Oa
    dCt = Ct * uCt - D * Ct

    dT = a * uCt * Ct / q[0] - D * T - Oa/q[2] * r[1] * T / (T + K[2])
    return [dR1, dR2, dT, dCt, dOa, uCt, uOa]

def rates(y,t,D):
    f = model(y,t,D)
    return(f[-2],f[-1])

def rates_chem(y,t,D):
    f = model(y,t,D)
    return(f[3],f[4])


def simulations():
    Ds = np.arange(0.0, 0.11, 0.01)
    for D in Ds:
        y = odeint(model, y0, xs, args=(D,))

        Oa = pd.DataFrame()
        Oa['N'], Oa['species'], Oa['x'] = y[:, 4], 'Oa', xs

        Ct = pd.DataFrame()
        Ct['N'], Ct['species'], Ct['x'] = y[:, 3], 'Ct', xs

        df = pd.concat([Ct, Oa])
        fig = px.line(df, x='x', y='N', color='species')
        fig.show()

def simulation(D):
    y = odeint(model, y0, xs, args=(D,))

    Oa = pd.DataFrame()
    Oa['N'], Oa['species'], Oa['x'] = y[:, 4], 'Oa', xs

    Ct = pd.DataFrame()
    Ct['N'], Ct['species'], Ct['x'] = y[:, 3], 'Ct', xs

    df = pd.concat([Ct, Oa])
    fig = px.line(df, x='x', y='N', color='species')
    fig.show()

def alle_effect():
    df = pd.DataFrame(columns=['N','u','species'])
    D = 0
    y = odeint(model, y0, xs, args=(D,))
    for j,i in enumerate(y):
        N = i[3] + i[4]
        uCt,uOa = rates(i,xs[j],D)
        df.loc[len(df)] = [N,uCt,'Ct']
        df.loc[len(df)] = [N,uOa,'Oa']

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

