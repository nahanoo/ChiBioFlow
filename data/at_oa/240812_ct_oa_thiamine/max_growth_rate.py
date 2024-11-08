import pandas as pd
from chibio_parser import fluorescence_paresr
import numpy as np
import plotly.express as px
from matplotlib import pyplot as plt
from scipy.stats import linregress


df = fluorescence_paresr("at_oa/240812_ct_oa_thiamine")
t1, t2 = 65, 80
df = df[(df["exp_time"] > t1) & (df["exp_time"] < t2)]


vs = []
# df = pd.read_csv("max_growth_rate.csv")
for r in ["M0", "M1", "M2"]:
    tmp = df[df["reactor"] == r]
    time = tmp["exp_time"].to_numpy()
    time -= time[0]
    data = tmp["od_measured"].to_numpy()
    slope, b = linregress(time, np.log(data))[:2]
    fit = [data[0] * np.exp(slope * t) for t in time]
    plt.plot(time, data, label=r)
    plt.plot(time, fit, label="Model" + r, linestyle="--")
    print(slope)
    vs.append(0.54 + slope)

plt.legend()
plt.show()
print(vs)
