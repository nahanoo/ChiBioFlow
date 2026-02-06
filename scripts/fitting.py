import pandas as pd
from data_parser import get_od_chemostats
import plotly.graph_objects as go
from style import *
import curveball
import numpy as np
from scipy.integrate import solve_ivp
import lmfit


def plot_od_chemostats(df):
    fig = go.Figure()
    for i, (_, df) in enumerate(df.groupby("name")):
        fig.add_trace(
            go.Scatter(
                x=df["time"][::5],
                y=df["OD"][::5],
                mode="markers",
                name=df["name"].iloc[0],
                marker=dict(
                    symbol="circle", size=6, color=list(colors_metabolites.values())[i]
                ),
                showlegend=False,
            )
        )
    fig = style_plot(fig)
    return fig


def plot_max_growth_rate(fig, x_means, rs, i):
    fig.add_trace(
        go.Scatter(
            x=x_means,
            y=rs,
            mode="markers",
            name=f"Replicate {i+1}",
            marker=dict(
                symbol="diamond", size=10, color=list(colors_metabolites.values())[i]
            ),
            showlegend=False,
        )
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]"), yaxis=dict(title="Growth rate [1/h]")
    )
    return fig


def plot_chemostat_fit(fig, t, N, i):
    fig.add_trace(
        go.Scatter(
            x=t,
            y=N,
            mode="lines",
            name="Simulated N",
            marker=dict(color=list(colors_metabolites.values())[i]),
            showlegend=False,
        )
    )
    fig.update_layout(xaxis=dict(title="Time [h]"), yaxis=dict(title="OD600"))
    return fig


def add_legend(fig, r, K, i):
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"r: {r:.3f} h⁻¹<br>Km: {K:.5f} mM",
            line=dict(color=list(colors_metabolites.values())[i]),
        )
    )
    return fig


def growth_rate_model(x, y, window=5, D=0.15):
    ln_y = np.log(y)
    rs = []
    x_means = []
    # Set time to 0 if clipped
    x = x - x[0]
    for i in range(len(x)):
        t0 = x[i]
        t1 = t0 + window
        mask = (x >= t0) & (x <= t1)
        if mask.sum() < 30:
            continue
        x_window = x[mask]
        ln_y_window = ln_y[mask]
        k = np.polyfit(x_window, ln_y_window, 1)[0]
        r = 0.15 + k
        rs.append(r)
        x_mean = x_window.mean()
        x_means.append(x_mean)
    return np.max(rs), (x_means, rs)


def chemostat_model(t, y, r, Km, Y, M, D):
    N, R = y
    dNdt = (r * R / (Km + R)) * N - D * N
    dRdt = D * (M - R) - (1 / Y) * (r * R / (Km + R)) * N
    return [dNdt, dRdt]


def simulate_model(t, r, Km, Y, M, D, N0, R0):
    sol = solve_ivp(
        fun=lambda t, y: chemostat_model(t, y, r, Km, Y, M, D),
        t_span=(t[0], t[-1]),
        y0=[N0, R0],
        t_eval=t,
        method="LSODA",
    )
    return sol.y


def fit_km(t, y, r, Km, Y, M, D, N0, R0):
    params = lmfit.Parameters()
    params.add("Km", value=Km, min=0.0001, max=20)

    def residual(p):
        Km = p["Km"].value
        N, _ = simulate_model(t, r, Km, Y, M, D, N0, R0)
        return N - y

    return lmfit.minimize(residual, params, method="least_squares")


def fit_max_growth_rate(df):
    fig = go.Figure()
    for i, (_, rep) in enumerate(df.groupby("name")):
        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()
        r, (x_means, r_windows) = growth_rate_model(x, y)
        plot_max_growth_rate(fig, x_means, r_windows, i)
    fig = style_plot(fig)
    fig.write_image("plots/fitting/oa_chemostat_growth_rates.svg")


def fit_oa(df):
    df = df[df["species"] == "Oa"]
    df = df[(df["time"] >= 1) & (df["time"] <= 60)]
    fig_N = plot_od_chemostats(df)
    fig_r = go.Figure()
    rs, Kms = [], []
    for i, (_, rep) in enumerate(df.groupby("name")):
        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()
        r, (x_means, r_windows) = growth_rate_model(x, y)
        r = np.mean(r_windows[: 60 * 25])
        # r = 0.19
        Y = max(y) / 7.5
        N0 = np.mean(y[:60])
        fit = fit_km(
            t=x,
            y=y,
            r=r,
            Km=0.02,
            Y=Y,
            M=7.5,
            D=0.15,
            N0=N0,
            R0=7.5,
        )
        Km_fitted = fit.params["Km"].value
        rs.append(r)
        Kms.append(Km_fitted)
        N, R = simulate_model(x, r, Km_fitted, Y, 7.5, 0.15, N0, 7.5)
        fig_N = plot_chemostat_fit(fig_N, x, N, i)
        fig_N = add_legend(fig_N, r, Km_fitted, i)
        fig_r = plot_max_growth_rate(fig_r, x_means, r_windows, i)

    fig_N = style_plot(fig_N, font_size=10, marker_size=1.5)
    fig_N.update_layout(width=300, height=150)
    fig_N.write_image("plots/fitting/oa_chemostat_fit.svg")
    fig_r = style_plot(fig_r, font_size=10, marker_size=1.5)
    fig_r.update_layout(width=150, height=150)
    fig_r.write_image("plots/fitting/oa_chemostat_growth_rates.svg")
    print(f"OA growth rate: {np.mean(rs):.3f} ± {np.std(rs):.3f}")
    print(f"OA Km: {np.mean(Kms):.3f} ± {np.std(Kms):.3f}")


def fit_ct(df):
    # df = pd.read_csv("../data/od_chemostats.csv")
    df = df[df["species"] == "Ct"]
    df = df[(df["time"] >= 1) & (df["time"] <= 12)]
    fig_N = plot_od_chemostats(df)
    fig_r = go.Figure()
    rs, Kms = [], []
    for i, (_, rep) in enumerate(df.groupby("name")):
        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()
        r, (x_means, r_windows) = growth_rate_model(x, y)
        r = r_windows[0]
        Y = (max(y) - 0.02) / 7.5
        N0 = np.mean(y[:50])
        fit = fit_km(
            t=x,
            y=y,
            r=r,
            Km=0.02,
            Y=Y,
            M=7.5,
            D=0.15,
            N0=N0,
            R0=7.5,
        )
        Km_fitted = fit.params["Km"].value
        N, R = simulate_model(x, r, Km_fitted, Y, 7.5, 0.15, N0, 7.5)
        fig_N = plot_chemostat_fit(fig_N, x, N, i)
        fig_N = add_legend(fig_N, r, Km_fitted, i)
        fig_r = plot_max_growth_rate(fig_r, x_means, r_windows, i)
        rs.append(r)
        Kms.append(Km_fitted)

    fig_N = style_plot(fig_N, font_size=10, marker_size=1.5)
    fig_N.update_layout(width=300, height=150)
    fig_N.write_image("plots/fitting/ct_chemostat_fit.svg")
    fig_r = style_plot(fig_r, font_size=10, marker_size=1.5)
    fig_r.update_layout(width=150, height=150)
    fig_r.write_image("plots/fitting/ct_chemostat_growth_rates.svg")
    print(f"CT growth rate: {np.mean(rs):.3f} ± {np.std(rs):.3f}")
    print(f"CT Km: {np.mean(Kms):.3f} ± {np.std(Kms):.3f}")


def normalized_growth_rates():
    fig = go.Figure()
    df = get_od_chemostats(write_excel=False)
    df = df[df["species"] == "Oa"]
    df = df[(df["time"] >= 1) & (df["time"] <= 60)]
    for i, (_, rep) in enumerate(df.groupby("name")):
        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()
        r, (x_means, r_windows) = growth_rate_model(x, y)
        fig.add_trace(
            go.Scatter(
                x=x_means,
                y=np.array(r_windows) / 0.18,
                mode="markers",
                marker=dict(color=colors["oa"]),
                showlegend=False,
            )
        )
    df = get_od_chemostats(write_excel=False)
    # df = pd.read_csv("../data/od_chemostats.csv")
    df = df[df["species"] == "Ct"]
    df = df[(df["time"] >= 1) & (df["time"] <= 12)]
    for i, (_, rep) in enumerate(df.groupby("name")):
        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()
        r, (x_means, r_windows) = growth_rate_model(x, y)
        fig.add_trace(
            go.Scatter(
                x=x_means,
                y=np.array(r_windows) / 0.45,
                mode="markers",
                marker=dict(color=colors["ct"]),
                showlegend=False,
            )
        )
    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="Normalized growth rate"),
        width=150,
        height=150,
    )

    fig = style_plot(fig, font_size=12, marker_size=2)
    fig.write_image("plots/fitting/chemostat_normalized_growth_rates.svg")


df = get_od_chemostats(write_excel=False)

fit_oa(df)
fit_ct(df)
normalized_growth_rates()
