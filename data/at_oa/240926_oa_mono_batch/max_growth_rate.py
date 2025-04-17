import pandas as pd
import numpy as np
from chibio_parser import fluorescence_paresr, cfu_parser
from os.path import join
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import linregress
from scipy.signal import savgol_filter
from scipy.integrate import odeint
import scipy.optimize as optimize
import plotly.graph_objects as go
from style_plot import *


def monod(y, t, v, Km, q):
    n, c = y
    dndt = v * c / (Km + c) * n
    dcdt = -v * c / (Km + c) * n / q
    return np.array([dndt, dcdt])


def simulate_monod(Km, v, t, q, n, c0, n0):
    y = odeint(monod, [n0, c0], t, args=(v, Km[0], q))
    return np.sum((n - y[:, 0])**2)


def get_Km(t, series, c0, n0, v, q):
    Km = optimize.minimize(
        simulate_monod,
        [0.1],
        args=(
            v,
            t,
            q,
            series,
            c0,
            n0,
        ),
        bounds=((0, 100), ),
    ).x
    return Km[0]


def plot_oa_chibio_batch(plot_log=False):
    dir = join("/", "home", "eric", "ChiBioFlow", "data")
    df = fluorescence_paresr(join(dir, "at_oa/240926_oa_mono_batch"))
    M0 = df[df["reactor"] == "M0"]
    M1 = df[df["reactor"] == "M1"]
    M2 = df[df["reactor"] == "M2"]
    t1 = 3.24
    t2 = 3.6
    for i, t in zip(M2.index, M2["exp_time"]):
        if t >= t1:
            it1 = i
            break
    for i, t in zip(M2.index, M2["exp_time"]):
        if t >= t2:
            it2 = i
            break
    raw1, raw2 = M1.loc[it1, "FP3_emit1"], M1.loc[it2, "FP3_emit1"]
    M1.loc[it2:, "FP3_emit1"] = M1.loc[it2:, "FP3_emit1"] - (raw2 - raw1)
    M1.loc[it1:it2, "FP3_emit1"] = raw1
    M1.loc[:, "exp_time"] = np.linspace(0, M1.iloc[-1]["exp_time"], len(M1))
    t1 = 10.5
    t2 = 11
    for i, t in zip(M2.index, M2["exp_time"]):
        if t >= t1:
            it1 = i
            break
    for i, t in zip(M2.index, M2["exp_time"]):
        if t >= t2:
            it2 = i
            break
    raw1, raw2 = M2.loc[it1, "FP3_emit1"], M2.loc[it2, "FP3_emit1"]
    M2.loc[it2:, "FP3_emit1"] = M2.loc[it2:, "FP3_emit1"] - (raw2 - raw1)
    M2.loc[it1:it2, "FP3_emit1"] = raw1
    M2.loc[:, "exp_time"] = np.linspace(0, M2.iloc[-1]["exp_time"], len(M2))
    reactors = [M0, M1, M2]
    od_standard = 0.45
    fig = make_subplots(rows=1, cols=1)
    legend = True
    y_windows, x_windows = [], []
    ys, xs = [], []
    exp_phase = len(M1[M1["exp_time"] < 8])
    stat_phase = len(M1[M1["exp_time"] < 14])
    for i, r in enumerate(reactors):
        x, y = r["exp_time"].to_numpy(), r["FP3_emit1"].to_numpy()
        x = x[:stat_phase]
        y = y[:len(x)]
        y = od_standard / y[-1] * y
        x_window = x[:exp_phase]
        y_window = y[:len(x_window)]
        xs.append(x)
        ys.append(y)
        window_size = 11
        # y = savgol_filter(y, window_length=window_size, polyorder=2)
        # x, y = x[::window_size], y[::window_size]
        if plot_log:
            plot_y = np.log(y)
        else:
            plot_y = y
        fig.add_trace(
            go.Scatter(
                x=x[::10],
                y=plot_y[::10],
                line=dict(color=colors['oa']),
                opacity=1,
                name="OD",
                showlegend=legend,
            ),
            row=1,
            col=1,
        )
        legend = False
        y_windows.append(y_window)
        x_windows.append(x_window)

    x_window_avg = np.average(np.array(x_windows), axis=0)
    y_window_avg = np.average(np.array(y_windows), axis=0)
    r, intercept = linregress(x_window_avg, np.log(y_window_avg))[:2]
    fit = [y[0] * np.exp(r * t) for t in x_window_avg]

    y1 = 0.07
    # fig.show()
    q = (0.46 - y1) / 7.5
    xs = np.average(np.array(xs), axis=0)
    ys = np.average(np.array(ys), axis=0)
    Km = get_Km(xs, ys, 7.5, ys[0], r, q)
    Y = odeint(monod, [y1, 7.5], xs, args=(r, Km, q))
    print('v', r, 'Km', Km)

    if plot_log:
        plot_y = np.log(Y[:, 0])
    else:
        plot_y = Y[:, 0]
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=plot_y,
            line=dict(color=colors['total'], dash="dot"),
            name="Fit",
            showlegend=True,
        ),
        row=1,
        col=1,
    )
    fig.update_xaxes(title="Time [h]")
    fig.update_yaxes(title="OD")
    return fig
