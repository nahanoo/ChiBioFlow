import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

f = 'citrate_fluorescence.csv'
fa = 'citrate_detection_test_diluted_absorbance.csv'
ft = 'standard_tests.csv'
df = pd.read_csv(ft)
standard_wells = ['F2', 'F3','F4', 'F5']
citrate_wells = ['B2','B3','D2','D3','D4']
standards = {w: None for w in standard_wells}
concentrations = [400, 240,120, 0]

for standard in standards:
    standards[standard] = np.average(df[standard]) - np.average(df['F5'])


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
        samples[sample] = (np.average(df[sample]) - np.average(df['F5'])) / m
    return samples


s = sample_conc(citrate_wells)
