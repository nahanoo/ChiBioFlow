import plotly.graph_objects as go
from scipy.integrate import odeint
from sympy import symbols, solve, diff, Matrix, Eq
from alle_effect import simulation
import numpy as np
from scipy.linalg import eigvals
import matplotlib.pyplot as plt
from os.path import join
w = 400
h = 250
dpi = 100
plt.rc('font', size=10) 


rCt = 0.24
rOa = 0.43
KC = 7
T_max = 0.0015
KT = T_max / 2
M = symbols('M')
a = 2 * T_max / M
qCt = 0.065
qOa = 0.11
qOaT = 1/T_max
KT = T_max / 2


Ct, Oa, C, T, D = symbols('Ct,Oa,C,T,D')

JCt = rCt * C / (C + KC)
JOa = rOa * C / (C + KC) * T / (T + KT)

dCt = JCt * Ct - D * Ct
dOa = JOa * Oa - D * Oa
dCom = JCt * Ct + JOa * Oa - D * (Ct + Oa)

dC = D * (M - C) - JCt * Ct / qCt - JOa * Oa / qOa
dT = a * JCt * Ct / qCt - JOa * Oa / qOaT - D * T


def growth_rates():
    C_values = np.linspace(0, 10, 50)
    D = 0.1
    Cf = solve(Eq(D, JCt))[0]
    Tf = solve(Eq(D, JOa), T)[0]
    Ts = [Tf.subs({'C': C}) for C in C_values]
    Cs = [Cf.subs({'C': C}) for C in C_values]
    T_values = np.linspace(0, 0.00148, 50)
    C_grid, T_grid = np.meshgrid(C_values, T_values)
    dJCt = rCt * C_grid / (C_grid + KC)
    dJOa = rOa * C_grid / (C_grid + KC) * T_grid / (T_grid + KT)


    plt.figure(figsize=(w*2/3/dpi,h/dpi),dpi=dpi)
    plt.contourf(C_values,T_values,dJCt,levels=100,cmap='RdYlBu')
    plt.colorbar(label='J [1/h]')
    line = plt.contour(C_values,T_values,dJCt,levels=[0.1],colors='black')
    label = ['$\mu$ = 0.1']
    plt.legend(line.collections,label,loc=1)
    plt.xlabel('Citric acid [mM]')
    plt.ylabel('Thiamine [mM]')
    plt.title(label='C. testosteroni',style='italic')
    plt.savefig(join('..','notebooks','mid_thesis_plots','jct.svg'), format='svg', bbox_inches='tight')


    plt.figure(figsize=(w*2/3/dpi,h/dpi),dpi=dpi)
    plt.contourf(C_values,T_values,dJOa,levels=100,cmap='RdYlBu')
    plt.colorbar(label='J [1/h]')
    line = plt.contour(C_values,T_values,dJOa,levels=[0.1],colors='black')
    label = ['$\mu$ = 0.1']
    plt.legend(line.collections,label,loc=1)
    plt.xlabel('Citric acid [mM]')
    plt.ylabel('Thiamine [mM]')
    plt.title(label='B. anthropi',style='italic')
    plt.savefig(join('..','notebooks','mid_thesis_plots','joa.svg'), format='svg', bbox_inches='tight')


def stable():
    D = 0.1
    Cf = solve(Eq(D, JCt))[0]
    Cs = np.linspace(0, 10, 50)
    Tf = solve(Eq(D, JOa), T)[0]
    Ts = []
    Cx = []
    for i,c in enumerate(Cs):
        t = Tf.subs({'C': c})
        if (t < 0.00148) & (t > 0):
            Ts.append(t)
            Cx.append(c)

    plt.figure(figsize=(w*2/3/dpi,h/dpi),dpi=dpi)
    plt.plot(Cx, Ts,color='#D95F02',label='B. anthropi')
    plt.plot([Cf]*len(Ts), Ts,color='#7570B3',label='C. testosteroni')
    plt.xlabel('Citric acid [mM]')
    plt.ylabel('Thiamine [mM]')
    plt.title(label='D = 0.1 [1/h]')
    legend = plt.legend()
    legend.get_texts()[0].set_fontstyle('italic')
    legend.get_texts()[1].set_fontstyle('italic')
    plt.savefig(join('..','notebooks','mid_thesis_plots','jctjoa.svg'), format='svg', bbox_inches='tight')

