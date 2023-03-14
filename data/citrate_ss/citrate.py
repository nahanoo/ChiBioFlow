import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

f = 'triplicates.csv'
df = pd.read_csv(f)
standard_wells = ['G1', 'G2', 'G3', 'G4']
citrate_wells = ['F'+str(i) for i in range(1, 13)]
standard_curve_samples = ['E1', 'E2', 'E3', 'E4']
standards = {w: None for w in standard_wells}
concentrations = [400, 240, 120, 0]

for standard in standards:
    fs = 'standard_curves.csv'
    dfs = pd.read_csv(fs)
    standards[standard] = np.average(dfs[standard]) - np.average(dfs['G4'])


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
        samples[sample] = (np.average(df[sample]) -
                           np.average(dfs['F4'])) / m
    return samples


s = sample_conc(citrate_wells)
