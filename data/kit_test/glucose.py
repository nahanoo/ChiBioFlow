import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

f = '../continuous_mono_test/glucose_steady_state_correct.csv'
df = pd.read_csv(f)
standard_wells = ['F1', 'F2','F3', 'F4','F5','F6']
glucose_stock = ['D1','D2','D3','D4']
samples = ['E1','E2','E3','E4']
standards = {w: None for w in standard_wells}
concentrations = [0.16,0.08,0.04,0.02,0.01,0]

for standard in standards:
    standards[standard] = np.average(df[standard]) - np.average(df['F6'])


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
        samples[sample] = (np.average(df[sample]) - np.average(df['G2'])) / m
    return samples


s = sample_conc(glucose_stock)
