import plotly.graph_objects as go
from style import *
from scipy.stats import linregress
from sympy import symbols, Eq, solve, Min
import numpy as np

cal = {0.01: 0.001, 0.025: 0.016, 0.05: 0.029, 0.1: 0.053}
media_raw = {"b1": 0.02, "b2": 0.047, "b3": 0.023}
ct_raw = {"ct1": -0.019, "ct2": -0.034, "ct3": -0.017}
oa_raw = {"oa1": -0.03, "oa2": -0.031, "oa3": -0.022}
x = list(cal.values())
y = list(cal.keys())
slope, intercept, r_value, p_value, std_err = linregress(x, y)
fit = [slope * i + intercept for i in np.linspace(0, x[-1], 20)]
fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y, mode="markers", name="Calibration"))
fig.add_trace(go.Scatter(x=np.linspace(0, x[-1], 20), y=fit, mode="lines", name="Fit"))
fig.update_layout(
    xaxis=dict(title="Difference in absorbance", range=[0, 0.065], dtick=0.02),
    yaxis=dict(title="Acetate mg/mL", range=[0, 0.12], dtick=0.02),
)
fig = style_plot(fig, marker_size=10)
fig.write_image("calibration.svg")


media = {
    key: (r * slope + intercept) * 4 / 60.05 * 1000 for key, r in media_raw.items()
}
ct = {key: (r * slope + intercept) / 60.05 * 1000 for key, r in ct_raw.items()}
oa = {key: (r * slope + intercept) / 60.05 * 1000 for key, r in oa_raw.items()}
"""measurements = {"c1": 1.313, "c2": 1.534, "c3": 1.16}
concentrations = dict()
for c, r in measurements.items():
    concentrations[c] = (r * slope + intercept) * 1e-6 / 0.266e-3 / 59.044 * 4 * 1e3

J = eq_J1_1
J = Eq(0.15, J.rhs)
K = solve(J, K1_1)[0]
print(K.subs({R1: concentrations["c2"], v1_1: 0.18}))"""
