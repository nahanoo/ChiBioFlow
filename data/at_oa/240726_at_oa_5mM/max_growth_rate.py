import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

df = pd.read_excel("at_oa_growth_curves.xlsx", sheet_name="Sheet2")
xs = np.linspace(0, 42.3, len(df))

at_well, oa_well = "C4", "F4"
at = df[["Time", at_well]]
oa = df[["Time", oa_well]]
r_at = np.gradient(at[at_well]) / at[at_well] * 3
r_oa = np.gradient(oa[oa_well]) / oa[oa_well] * 3
max_r_at = max(r_at)
max_r_oa = max(r_oa)


def plot_growth_curve():
    plt.plot(xs, at[at_well], label="At")
    plt.plot(xs, oa[oa_well], label="Oa")
    plt.legend()
    plt.show()


plot_growth_curve()
