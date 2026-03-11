import pandas as pd
from data_parser import get_od_chemostats
import plotly.graph_objects as go
from style import *
import curveball
import numpy as np
from scipy.integrate import solve_ivp
import lmfit
from plotly.subplots import make_subplots


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


def fit_km(t, y, r, Km, Y, M, D, N0, R0, fit_start=30):
    params = lmfit.Parameters()
    params.add("Km", value=Km, min=0.0001, max=20)

    fit_mask = t >= fit_start

    def residual(p):
        Km = p["Km"].value
        N, _ = simulate_model(t, r, Km, Y, M, D, N0, R0)
        return N[fit_mask] - y[fit_mask]

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
    rs, Kms, Ys = [], [], []
    for i, (_, rep) in enumerate(df.groupby("name")):
        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()
        r, (x_means, r_windows) = growth_rate_model(x, y)
        r = np.mean(r_windows[: 60 * 25])
        # r = 0.19
        N0 = np.mean(y[:60])
        Y = max(y - 0.02) / 7.5
        Ys.append(Y)
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
            fit_start=30,
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
    print(f"OA yield: {np.mean(Ys):.3f} ± {np.std(Ys):.3f}")


def fit_ct(df):
    # df = pd.read_csv("../data/od_chemostats.csv")
    df = df[df["species"] == "Ct"]
    df = df[(df["time"] >= 1) & (df["time"] <= 12)]
    fig_N = plot_od_chemostats(df)
    fig_r = go.Figure()
    rs, Kms, Ys = [], [], []
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
            fit_start=5,
        )
        Km_fitted = fit.params["Km"].value
        N, R = simulate_model(x, r, Km_fitted, Y, 7.5, 0.15, N0, 7.5)
        fig_N = plot_chemostat_fit(fig_N, x, N, i)
        fig_N = add_legend(fig_N, r, Km_fitted, i)
        fig_r = plot_max_growth_rate(fig_r, x_means, r_windows, i)
        rs.append(r)
        Kms.append(Km_fitted)
        Ys.append(Y)
    fig_N = style_plot(fig_N, font_size=10, marker_size=1.5)
    fig_N.update_layout(width=300, height=150)
    fig_N.write_image("plots/fitting/ct_chemostat_fit.svg")
    fig_r = style_plot(fig_r, font_size=10, marker_size=1.5)
    fig_r.update_layout(width=150, height=150)
    fig_r.write_image("plots/fitting/ct_chemostat_growth_rates.svg")
    print(f"CT growth rate: {np.mean(rs):.3f} ± {np.std(rs):.3f}")
    print(f"CT Km: {np.mean(Kms):.3f} ± {np.std(Kms):.3f}")
    print(f"CT yield: {np.mean(Ys):.3f} ± {np.std(Ys):.3f}")


def parse_params():
    df = dict(pd.read_csv("parameters.csv"))
    params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
    p = params
    return p


def thiamine_supply_rhs(t, y, p):
    Ct, Oa, R, T = y

    JCt = p["v1_1"] * R / (R + p["K1_1"])
    JOa = p["v2_1"] * R / (R + p["K2_1"]) * T / (T + p["K2_3"])

    dCt = JCt * Ct - p["D"] * Ct
    dOa = JOa * Oa - p["D"] * Oa
    dR = -JCt * Ct / p["q1_1"] - JOa * Oa / p["q2_1"] + p["D"] * p["M1"] - p["D"] * R
    dT = -JOa * Oa / p["q2_3"] - p["D"] * T + p["M3"] * p["D"]

    return [dCt, dOa, dR, dT]


def simulate_thiamine_supply(x, p, y0):
    x = np.asarray(x, dtype=float)
    sol = solve_ivp(
        fun=lambda t, y: thiamine_supply_rhs(t, y, p),
        t_span=(x[0], x[-1]),
        y0=y0,
        t_eval=x,
        method="LSODA",
    )
    return sol.y


def fit_km_thiamine_1nM(
    x,
    y,
    p,
    y0,
    K2_3_init=1.0,
    K2_3_min=1e-4,
    K2_3_max=1e3,
    use_log=False,
):
    """
    Fit p['K2_3'] from a single 1 nM thiamine Oa trajectory.

    Units:
    - T in nM
    - K2_3 in nM
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    if use_log:
        mask &= y > 0

    x = x[mask]
    y = y[mask]

    params = lmfit.Parameters()
    params.add("K2_3", value=K2_3_init, min=K2_3_min, max=K2_3_max)

    def residual(pars):
        p_fit = p.copy()
        p_fit["K2_3"] = pars["K2_3"].value

        sim = simulate_thiamine_supply(x, p_fit, y0)
        y_sim = sim[1]  # Oa

        if use_log:
            y_sim = np.clip(y_sim, 1e-12, None)
            return np.log(y_sim) - np.log(y)

        return y_sim - y

    fit = lmfit.minimize(residual, params, method="least_squares")
    K2_3_fitted = fit.params["K2_3"].value
    return K2_3_fitted, fit


def plot_thiamine_fit_1nM(x, y, p, y0, K2_3_fitted, outname):
    p_fit = p.copy()
    p_fit["K2_3"] = K2_3_fitted

    sim = simulate_thiamine_supply(x, p_fit, y0)
    y_sim = sim[1]  # Oa

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(size=6, color=colors["1 nM thiamine"]),
            name="1 nM thiamine data",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y_sim,
            mode="lines",
            line=dict(color=colors["1 nM thiamine"]),
            name=f"fit K2_3={K2_3_fitted:.3g} nM",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD600"),
        width=300,
        height=200,
    )
    fig = style_plot(fig, font_size=10, marker_size=3, line_thickness=2)
    fig.write_image(outname)


def fit_thiamine():
    colors = {
        "0 nM thiamine": "#1f77b4",
        "0.01 nM thiamine": "#ff7f0e",
        "0.1 nM thiamine": "#2ca02c",
        "1 nM thiamine": "#d62728",
        "10 nM thiamine": "#9467bd",
        "100 nM thiamine": "#8c564b",
        "1000 nM thiamine": "#e377c2",
        "10000 nM thiamine": "#7f7f7f",
        "ct": "#7570B3",
        "oa": "#D95F02",
    }

    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/metadata.csv"
    )
    df = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_oa_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
        & (df["comments"] == "1 nM thiamine")
    ].copy()

    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/measurements.csv"
    )

    p = parse_params()
    p["D"] = 0

    for _, row in df.iterrows():
        lg = row["linegroup"]

        x = data[f"{lg}_time"].to_numpy()
        y = data[f"{lg}_measurement"].to_numpy()

        mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
        x = x[mask]
        y = y[mask]

        # Oa monoculture at 1 nM thiamine
        y0 = [0.0, 0.01, 5.0, 1.0]  # [Ct0, Oa0, R0, T0]

        K2_3_fitted, fit = fit_km_thiamine_1nM(
            x=x,
            y=y,
            p=p,
            y0=y0,
            K2_3_init=1.0,  # nM
            K2_3_min=1e-4,  # nM
            K2_3_max=1e3,  # nM
            use_log=False,
        )

        print("Fitted K2_3 [nM]:", K2_3_fitted)
        print(lmfit.fit_report(fit))

        plot_thiamine_fit_1nM(
            x=x,
            y=y,
            p=p,
            y0=y0,
            K2_3_fitted=K2_3_fitted,
            outname="plots/fitting/oa_thiamine_fit_1nM.svg",
        )


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


def plot_od_and_normalized_growth(df, species_name, t_min, t_max, outname):
    df = df[df["species"] == species_name].copy()
    df = df[(df["time"] >= t_min) & (df["time"] <= t_max)]

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    for i, (_, rep) in enumerate(df.groupby("name")):
        rep = rep.sort_values("time")

        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()

        # avoid log problems in growth_rate_model
        mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
        x = x[mask]
        y = y[mask]

        if len(x) == 0:
            continue

        # OD trace
        fig.add_trace(
            go.Scatter(
                x=x[::5],
                y=y[::5],
                mode="markers",
                marker=dict(
                    symbol="circle",
                    size=6,
                    color=list(colors_metabolites.values())[i],
                ),
                name=rep["name"].iloc[0],
                opacity=0.5,
                showlegend=False,
            ),
            secondary_y=False,
        )

        # normalized growth rate
        r_max, (x_means, r_windows) = growth_rate_model(x, y)
        mu_norm = np.array(r_windows) / r_max

        fig.add_trace(
            go.Scatter(
                x=x_means[::30],
                y=mu_norm[::30],
                mode="lines",
                marker=dict(
                    symbol="diamond",
                    size=7,
                    color=list(colors_metabolites.values())[i],
                ),
                line=dict(dash="dot"),
                showlegend=False,
            ),
            secondary_y=True,
        )

    # dummy legend entries
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(symbol="circle", size=6, color="gray"),
            name="Experimental OD",
        ),
        secondary_y=False,
    )

    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                symbol="diamond",
                size=7,
                color=colors[species_name.lower()],
            ),
            name="Normalized growth rate (μ/rmax)",
        ),
        secondary_y=True,
    )
    fig.update_xaxes(title_text="Time [h]")
    fig.update_yaxes(title_text="OD600", secondary_y=False)
    fig.update_yaxes(title_text="Normalized growth rate", secondary_y=True)

    fig.update_layout(
        width=300,
        height=180,
        showlegend=False,
    )

    fig = style_plot(fig, font_size=10, marker_size=3, line_thickness=2)
    fig.write_image(outname)


def plot_net_growth_rate(df, species_name, t_min, t_max, outname, D=0.15):
    """
    Plot dlnN/dt = mu - D versus time for one species.
    This goes to ~0 at steady state, so it is easier to interpret than mu/rmax.
    """
    df = df[df["species"] == species_name].copy()
    df = df[(df["time"] >= t_min) & (df["time"] <= t_max)]

    fig = go.Figure()

    for i, (_, rep) in enumerate(df.groupby("name")):
        rep = rep.sort_values("time")

        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()

        mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
        x = x[mask]
        y = y[mask]

        if len(x) == 0:
            continue

        _, (x_means, mus) = growth_rate_model(x, y, D=D)
        net_growth = np.array(mus) - D  # dlnN/dt

        fig.add_trace(
            go.Scatter(
                x=x_means,
                y=net_growth,
                mode="lines+markers",
                line=dict(color=list(colors_metabolites.values())[i]),
                marker=dict(
                    symbol="circle",
                    size=5,
                    color=list(colors_metabolites.values())[i],
                ),
                name=rep["name"].iloc[0],
                showlegend=False,
            )
        )

    fig.add_hline(y=0, line=dict(color="gray", dash="dash"))

    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="Net growth rate dln(N)/dt [1/h]"),
        width=300,
        height=180,
        showlegend=False,
    )

    fig = style_plot(fig, font_size=10, marker_size=3, line_thickness=2)
    fig.write_image(outname)


def plot_net_growth_vs_od(df, species_name, t_min, t_max, outname, D=0.15):
    """
    Plot dlnN/dt = mu - D versus OD for one species.
    This is more informative for Km because OD is a proxy for cumulative substrate use
    if yield is similar between species.
    """
    df = df[df["species"] == species_name].copy()
    df = df[(df["time"] >= t_min) & (df["time"] <= t_max)]

    fig = go.Figure()

    for i, (_, rep) in enumerate(df.groupby("name")):
        rep = rep.sort_values("time")

        x = rep["time"].to_numpy()
        y = rep["OD"].to_numpy()

        mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
        x = x[mask]
        y = y[mask]

        if len(x) == 0:
            continue

        _, (x_means, mus) = growth_rate_model(x, y, D=D)
        net_growth = np.array(mus) - D  # dlnN/dt

        # interpolate OD at the window midpoints
        od_means = np.interp(x_means, x, y)

        fig.add_trace(
            go.Scatter(
                x=od_means,
                y=net_growth,
                mode="lines+markers",
                line=dict(color=list(colors_metabolites.values())[i]),
                marker=dict(
                    symbol="circle",
                    size=5,
                    color=list(colors_metabolites.values())[i],
                ),
                name=rep["name"].iloc[0],
                showlegend=False,
            )
        )

    fig.add_hline(y=0, line=dict(color="gray", dash="dash"))

    fig.update_layout(
        xaxis=dict(title="OD600"),
        yaxis=dict(title="Net growth rate dln(N)/dt [1/h]"),
        width=120,
        height=150,
        showlegend=False,
    )

    fig = style_plot(fig, font_size=10, marker_size=3, line_thickness=5)
    fig.write_image(outname)


def plot_net_growth_both(df):
    plot_net_growth_rate(
        df=df,
        species_name="Oa",
        t_min=1,
        t_max=60,
        outname="plots/fitting/oa_net_growth_rate.svg",
    )
    plot_net_growth_rate(
        df=df,
        species_name="Ct",
        t_min=1,
        t_max=12,
        outname="plots/fitting/ct_net_growth_rate.svg",
    )

    plot_net_growth_vs_od(
        df=df,
        species_name="Oa",
        t_min=1,
        t_max=60,
        outname="plots/fitting/oa_net_growth_vs_od.svg",
    )
    plot_net_growth_vs_od(
        df=df,
        species_name="Ct",
        t_min=1,
        t_max=12,
        outname="plots/fitting/ct_net_growth_vs_od.svg",
    )


"""df = get_od_chemostats(write_excel=False)

fit_oa(df)
fit_ct(df)
"""
