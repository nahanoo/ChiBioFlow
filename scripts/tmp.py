import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.integrate import solve_ivp
import lmfit

from data_parser import get_od_chemostats
from style import *


# ----------------------------
# Settings
# ----------------------------
D = 0.15
M = 7.5
Y_SHARED = None  # if None, estimate from each replicate; otherwise set a fixed value

OA_TIME_MIN = 1
OA_TIME_MAX = 60
CT_TIME_MIN = 1
CT_TIME_MAX = 12

OA_FIT_START = 30
CT_FIT_START = 5

GROWTH_WINDOW = 5
MIN_POINTS_PER_WINDOW = 10

OUTDIR = "plots/fitting"
os.makedirs(OUTDIR, exist_ok=True)


# ----------------------------
# Utilities
# ----------------------------
def _rep_color(i):
    return list(colors_metabolites.values())[i % len(colors_metabolites)]


def clean_xy(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
    return x[mask], y[mask]


def estimate_yield(y, M=M):
    # identical yield assumption can be enforced by replacing this with a constant
    # or by pooling both species and taking a single mean estimate
    return np.max(y) / M


# ----------------------------
# Plotting
# ----------------------------
def plot_od_condition(df, species_name, outname):
    fig = go.Figure()
    for i, (_, rep) in enumerate(df.groupby("name")):
        rep = rep.sort_values("time")
        fig.add_trace(
            go.Scatter(
                x=rep["time"][::5],
                y=rep["OD"][::5],
                mode="markers",
                name=rep["name"].iloc[0],
                marker=dict(
                    symbol="circle",
                    size=6,
                    color=_rep_color(i),
                ),
                showlegend=False,
            )
        )
    fig.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD600"),
    )
    fig = style_plot(fig, font_size=10, marker_size=1.5)
    fig.write_image(os.path.join(OUTDIR, outname))
    return fig


def add_growth_trace(fig, x_means, rs, color, name=None):
    fig.add_trace(
        go.Scatter(
            x=x_means,
            y=rs,
            mode="markers",
            marker=dict(symbol="diamond", size=7, color=color),
            name=name,
            showlegend=name is not None,
        )
    )
    return fig


def add_norm_growth_trace(fig, x, y, color, name=None):
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(symbol="circle", size=6, color=color),
            name=name,
            showlegend=name is not None,
        )
    )
    return fig


def add_fit_overlay(fig, t_fit, y_fit, y_sim_fit, color, label):
    fig.add_trace(
        go.Scatter(
            x=t_fit,
            y=y_fit,
            mode="markers",
            marker=dict(size=6, color=color, symbol="circle-open"),
            name=f"{label} data (fit region)",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t_fit,
            y=y_sim_fit,
            mode="lines",
            line=dict(color=color),
            name=f"{label} model",
            showlegend=False,
        )
    )
    return fig


# ----------------------------
# Growth rate estimation
# ----------------------------
def growth_rate_model(
    x, y, window=GROWTH_WINDOW, D=D, min_points=MIN_POINTS_PER_WINDOW
):
    x, y = clean_xy(x, y)

    ln_y = np.log(y)
    x0 = x[0]
    x_shifted = x - x0

    x_means = []
    mus = []

    for i in range(len(x_shifted)):
        t0 = x_shifted[i]
        t1 = t0 + window
        mask = (x_shifted >= t0) & (x_shifted <= t1)

        if mask.sum() < min_points:
            continue

        x_window = x_shifted[mask]
        ln_y_window = ln_y[mask]

        slope = np.polyfit(x_window, ln_y_window, 1)[0]
        mu = D + slope

        x_means.append(x_window.mean() + x0)
        mus.append(mu)

    x_means = np.asarray(x_means)
    mus = np.asarray(mus)

    if len(mus) == 0:
        raise ValueError("No valid growth-rate windows found.")

    return np.max(mus), (x_means, mus)


def choose_r_for_species(species_name, x_means, mus):
    x_means = np.asarray(x_means)
    mus = np.asarray(mus)

    if species_name == "Oa":
        # use early region before the kink
        mask = x_means < 20
        if mask.sum() > 0:
            return np.max(mus[mask])
        return np.max(mus)

    if species_name == "Ct":
        # Ct transition is early, so just use the max early mu
        mask = x_means < 5
        if mask.sum() > 0:
            return np.max(mus[mask])
        return np.max(mus)

    return np.max(mus)


# ----------------------------
# Chemostat model
# ----------------------------
def chemostat_model(t, y, r, Km, Y, M, D):
    N, R = y
    mu = r * R / (Km + R)
    dNdt = mu * N - D * N
    dRdt = D * (M - R) - (1 / Y) * mu * N
    return [dNdt, dRdt]


def simulate_model(t, r, Km, Y, M, D, N0, R0):
    sol = solve_ivp(
        fun=lambda t_, y_: chemostat_model(t_, y_, r, Km, Y, M, D),
        t_span=(t[0], t[-1]),
        y0=[N0, R0],
        t_eval=t,
        method="LSODA",
    )
    return sol.y


def fit_km(t, y, r, Km0, Y, M, D, N0, R0, fit_start):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    fit_mask = np.isfinite(t) & np.isfinite(y) & (y > 0) & (t >= fit_start)

    params = lmfit.Parameters()
    params.add("Km", value=Km0, min=1e-5, max=20)

    def residual(p):
        Km = p["Km"].value
        N_sim, _ = simulate_model(t, r, Km, Y, M, D, N0, R0)
        return N_sim[fit_mask] - y[fit_mask]

    fit = lmfit.minimize(residual, params, method="least_squares")
    return fit, fit_mask


# ----------------------------
# Additional proof metric
# ----------------------------
def od_at_norm_growth_threshold(x, y, threshold=0.8):
    x = np.asarray(x)
    y = np.asarray(y)

    below = np.where(y < threshold)[0]
    if len(below) == 0:
        return np.nan
    return x[below[0]]


# ----------------------------
# Per-species analysis
# ----------------------------
def analyze_species(df, species_name, tmin, tmax, fit_start, out_prefix):
    df = df[df["species"] == species_name].copy()
    df = df[(df["time"] >= tmin) & (df["time"] <= tmax)].copy()

    plot_od_condition(df, species_name, f"{out_prefix}_od.svg")

    fig_growth = go.Figure()
    fig_fit = go.Figure()

    rs = []
    kms = []
    threshold_ods = []

    all_norm_growth = []

    for i, (_, rep) in enumerate(df.groupby("name")):
        rep = rep.sort_values("time")

        x = rep["time"].to_numpy(dtype=float)
        y = rep["OD"].to_numpy(dtype=float)
        x, y = clean_xy(x, y)

        r_max, (x_means, mus) = growth_rate_model(x, y)
        r = choose_r_for_species(species_name, x_means, mus)

        Y = Y_SHARED if Y_SHARED is not None else estimate_yield(y, M=M)
        N0 = np.mean(y[: min(20, len(y))])
        R0 = M

        fit, fit_mask = fit_km(
            t=x,
            y=y,
            r=r,
            Km0=0.02,
            Y=Y,
            M=M,
            D=D,
            N0=N0,
            R0=R0,
            fit_start=fit_start,
        )

        Km_fitted = fit.params["Km"].value
        N_sim, R_sim = simulate_model(x, r, Km_fitted, Y, M, D, N0, R0)

        # growth-rate plot
        add_growth_trace(
            fig_growth,
            x_means,
            mus,
            color=_rep_color(i),
            name=None,
        )

        # fit plot: show full OD faintly, but overlay only the fit region
        fig_fit.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="markers",
                marker=dict(size=4, color=_rep_color(i)),
                opacity=0.25,
                showlegend=False,
                name=rep["name"].iloc[0],
            )
        )
        add_fit_overlay(
            fig_fit,
            t_fit=x[fit_mask],
            y_fit=y[fit_mask],
            y_sim_fit=N_sim[fit_mask],
            color=_rep_color(i),
            label=rep["name"].iloc[0],
        )

        # normalized growth rate
        mu_norm = mus / r
        od_interp = np.interp(x_means, x, y)
        all_norm_growth.append(
            pd.DataFrame(
                {
                    "species": species_name,
                    "name": rep["name"].iloc[0],
                    "time": x_means,
                    "OD": od_interp,
                    "mu": mus,
                    "mu_norm": mu_norm,
                }
            )
        )

        threshold_od = od_at_norm_growth_threshold(od_interp, mu_norm, threshold=0.8)

        rs.append(r)
        kms.append(Km_fitted)
        threshold_ods.append(threshold_od)

        print(
            f"{species_name} | {rep['name'].iloc[0]} | "
            f"r = {r:.3f} 1/h | Km = {Km_fitted:.5f} mM | "
            f"OD@mu/r<0.8 = {threshold_od:.4f} | chi2 = {fit.chisqr:.5g}"
        )

    # save growth-rate figure
    fig_growth.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="Growth rate [1/h]"),
        width=300,
        height=150,
    )
    fig_growth = style_plot(fig_growth, font_size=10, marker_size=1.5)
    fig_growth.write_image(os.path.join(OUTDIR, f"{out_prefix}_growth_rates.svg"))

    # save fit figure
    fig_fit.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="OD600"),
        width=300,
        height=150,
    )
    fig_fit = style_plot(fig_fit, font_size=10, marker_size=1.5)
    fig_fit.write_image(os.path.join(OUTDIR, f"{out_prefix}_km_fit_kink_only.svg"))

    norm_df = pd.concat(all_norm_growth, ignore_index=True)

    summary = {
        "species": species_name,
        "r_mean": np.mean(rs),
        "r_sd": np.std(rs),
        "Km_mean": np.mean(kms),
        "Km_sd": np.std(kms),
        "OD_threshold_mean": np.nanmean(threshold_ods),
        "OD_threshold_sd": np.nanstd(threshold_ods),
        "norm_df": norm_df,
    }

    print(f"{species_name}: r = {summary['r_mean']:.3f} ± {summary['r_sd']:.3f} 1/h")
    print(f"{species_name}: Km = {summary['Km_mean']:.5f} ± {summary['Km_sd']:.5f} mM")
    print(
        f"{species_name}: OD at mu/r < 0.8 = "
        f"{summary['OD_threshold_mean']:.4f} ± {summary['OD_threshold_sd']:.4f}"
    )

    return summary


# ----------------------------
# Combined normalized growth plot
# ----------------------------
def plot_normalized_growth_combined(norm_oa, norm_ct):
    fig_time = go.Figure()
    fig_od = go.Figure()

    for _, rep in norm_oa.groupby("name"):
        add_norm_growth_trace(
            fig_time,
            rep["time"],
            rep["mu_norm"],
            color=colors["oa"],
            name=None,
        )
        add_norm_growth_trace(
            fig_od,
            rep["OD"],
            rep["mu_norm"],
            color=colors["oa"],
            name=None,
        )

    for _, rep in norm_ct.groupby("name"):
        add_norm_growth_trace(
            fig_time,
            rep["time"],
            rep["mu_norm"],
            color=colors["ct"],
            name=None,
        )
        add_norm_growth_trace(
            fig_od,
            rep["OD"],
            rep["mu_norm"],
            color=colors["ct"],
            name=None,
        )

    # dummy legend
    fig_time.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(color=colors["oa"], size=7),
            name="Oa",
        )
    )
    fig_time.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(color=colors["ct"], size=7),
            name="Ct",
        )
    )

    fig_od.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(color=colors["oa"], size=7),
            name="Oa",
        )
    )
    fig_od.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(color=colors["ct"], size=7),
            name="Ct",
        )
    )

    fig_time.update_layout(
        xaxis=dict(title="Time [h]"),
        yaxis=dict(title="Normalized growth rate (μ/r)"),
        width=180,
        height=150,
    )
    fig_time = style_plot(fig_time, font_size=10, marker_size=1.5)
    fig_time.write_image(os.path.join(OUTDIR, "normalized_growth_rate_vs_time.svg"))

    fig_od.update_layout(
        xaxis=dict(title="OD600"),
        yaxis=dict(title="Normalized growth rate (μ/r)"),
        width=180,
        height=150,
    )
    fig_od = style_plot(fig_od, font_size=10, marker_size=1.5)
    fig_od.write_image(os.path.join(OUTDIR, "normalized_growth_rate_vs_od.svg"))


# ----------------------------
# Summary comparison plot
# ----------------------------
def plot_threshold_summary(summary_oa, summary_ct):
    fig = go.Figure()

    xs = ["Oa", "Ct"]
    ys = [
        summary_oa["OD_threshold_mean"],
        summary_ct["OD_threshold_mean"],
    ]
    errs = [
        summary_oa["OD_threshold_sd"],
        summary_ct["OD_threshold_sd"],
    ]
    cols = [colors["oa"], colors["ct"]]

    for x, y, err, c in zip(xs, ys, errs, cols):
        fig.add_trace(
            go.Scatter(
                x=[x],
                y=[y],
                mode="markers",
                marker=dict(size=10, color=c),
                error_y=dict(type="data", array=[err], visible=True),
                showlegend=False,
            )
        )

    fig.update_layout(
        xaxis=dict(title="Species"),
        yaxis=dict(title="OD at μ/r < 0.8"),
        width=180,
        height=150,
    )
    fig = style_plot(fig, font_size=10, marker_size=1.5)
    fig.write_image(os.path.join(OUTDIR, "od_at_norm_growth_drop_summary.svg"))


# ----------------------------
# Main
# ----------------------------
def main():
    df = get_od_chemostats(write_excel=False)

    oa = analyze_species(
        df=df,
        species_name="Oa",
        tmin=OA_TIME_MIN,
        tmax=OA_TIME_MAX,
        fit_start=OA_FIT_START,
        out_prefix="oa",
    )

    ct = analyze_species(
        df=df,
        species_name="Ct",
        tmin=CT_TIME_MIN,
        tmax=CT_TIME_MAX,
        fit_start=CT_FIT_START,
        out_prefix="ct",
    )

    plot_normalized_growth_combined(oa["norm_df"], ct["norm_df"])
    plot_threshold_summary(oa, ct)

    print("\nComparison")
    print(
        f"Km_Oa = {oa['Km_mean']:.5f} ± {oa['Km_sd']:.5f} mM | "
        f"Km_Ct = {ct['Km_mean']:.5f} ± {ct['Km_sd']:.5f} mM"
    )
    print(
        f"OD at μ/r < 0.8: Oa = {oa['OD_threshold_mean']:.4f} ± {oa['OD_threshold_sd']:.4f}, "
        f"Ct = {ct['OD_threshold_mean']:.4f} ± {ct['OD_threshold_sd']:.4f}"
    )
    print(
        "Interpretation: if Oa reaches a higher OD before normalized growth drops, "
        "that supports lower effective Km for Oa."
    )


if __name__ == "__main__":
    main()
