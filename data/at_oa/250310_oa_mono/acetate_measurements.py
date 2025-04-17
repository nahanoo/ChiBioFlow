import plotly.graph_objects as go
from style import *
from scipy.stats import linregress
from sympy import symbols, Eq, solve, Min
from equations import *


cal = {0.1: 1.507, 0.5: 1.122, 1: 1.007}
x = list(cal.values())
y = list(cal.keys())
slope, intercept, r_value, p_value, std_err = linregress(x, y)
fit = [slope * i + intercept for i in x]
fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y, mode="markers", name="Calibration"))
fig.add_trace(go.Scatter(x=x, y=fit, mode="lines", name="Fit"))
fig = style_plot(fig, marker_size=10)
fig.write_image("calibration.svg")

measurements = {"c1": 1.313, "c2": 1.534, "c3": 1.16}
concentrations = dict()
for c, r in measurements.items():
    concentrations[c] = (r * slope + intercept) * 1e-6 / 0.266e-3 / 59.044 * 4 * 1e3

J = eq_J1_1
J = Eq(0.15, J.rhs)
K = solve(J, K1_1)[0]
print(K.subs({R1: concentrations["c2"], v1_1: 0.18}))
