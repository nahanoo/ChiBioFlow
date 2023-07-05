import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
from os.path import join
from matplotlib import pyplot as plt
from scipy.integrate import odeint
import plotly.express as px


f = 'ss_30_3.csv'
df = pd.read_csv(f)
standard_wells = ['H1', 'H2', 'H3', 'H4']
citrate_wells = ['F'+str(i) for i in range(1, 13)]
standards = {w: None for w in standard_wells}
concentrations = [400, 240, 120, 0]

for standard in standards:
    standards[standard] = np.average(
        df[standard]) - np.average(df[standard_wells[-1]])


def linear_f(x, m, t):
    return t + m * x


def regression(plot=False):
    (m, t), cv = curve_fit(linear_f, concentrations, list(standards.values()))
    ys = [linear_f(x, m, t) for x in concentrations]
    if plot:
        plt.plot(concentrations, standards.values())
        plt.plot(concentrations, ys)
        plt.show()
    return m


def sample_conc(wells):
    m = regression()
    samples = {w: None for w in wells}
    for sample in samples:
        samples[sample] = np.average(df[sample]) / m
    return samples


def average(s):
    # Only works for triplicates
    for i in np.arange(0, 12, 3):
        av = np.average(list(s.values())[i:i+3])
        print(av)


def max_u():
    pass


f_c = pd.read_csv('ss_17_3_ct.csv')[['B3']]
f_o = pd.read_csv('oa_03_19_curves.csv')[['Citric acid']]
u_c = np.gradient(f_c['B3'], 1/6) / f_c['B3']
u_o = np.gradient(f_o['Citric acid'], 1/6) / f_o['Citric acid']
print('Ct', max(u_c))
print('Oa', max(u_o))
# return max(u_c), max(u_o)


def get_Y_cfus():
    ct_cfus = 2700000000.0
    oa_cfus = 13333333333.333334
    s_ct = 137.4976471744328 * 50.0
    f_ct = 196.88529428319006 * 50.0
    s_oa = 140.54647070640806 * 50.0
    f_oa = 202.25294134652674 * 50.0
    y_ct = ct_cfus / (f_ct - s_ct)
    y_oa = oa_cfus / (f_oa - s_oa)
    print('Ct', y_ct, 'Oa', y_oa)
    return y_ct, y_oa


def get_Y_OD():
    ct_cfus = 0.26
    oa_cfus = 0.17
    s_ct = 97.77705890574124 * 50.0
    f_ct = 159.4835295458599 * 50.0
    s_oa = 94.728235373766 * 50.0
    f_oa = 154.11588248252323 * 50.0
    y_ct = ct_cfus / (f_ct - s_ct)
    y_oa = oa_cfus / (f_oa - s_oa)
    print('Ct', y_ct, 'Oa', y_oa)
    return y_ct, y_oa


def get_K():
    u_c, u_o = 0.1, 0.1
    K_c = 137.4976471744328 * 50.0 * (u_c - 0.05) / 0.05
    K_o = 140.54647070640806 * 50.0 * (u_o - 0.05) / 0.05
    print('Ct', K_c, 'Oa', K_o)
    return K_c, K_o


M = 7974.1764772929955


def model(y, t, u, K, q, D):
    N = y[0]
    R = y[1]
    dN = u * R / (K + R) * N - D * N
    dR = D * M - D * R - N / \
        q * u * R / (K + R)
    return [dN, dR]


D = 0
u = 0.1
K = get_K()[0]
q = get_Y_OD()[0]
x = np.arange(0, 433/10, 0.1)
y = odeint(model, [0.1, M], x, args=(u, K, q, D))
OD = y[:, 0]
plt.plot(x, OD)
plt.plot(x, f_c['B3'])
# plt.show()


#s = sample_conc(citrate_wells)
