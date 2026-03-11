import pandas as pd
from data_parser import df, get_cfus, get_od_chemostats, ct_oa_plate_reader
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from chibio_parser import *
from models import *
from scipy.stats import linregress
from style import *
import scipy.stats as stats
import statsmodels.formula.api as smf
import curveball


def growth_curves_ct_oa():
    model = False
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250328_ct_oa_thiamine_gradient/data/metadata.csv"
    )
    df_ct = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
        & (df["species"] == "Comamonas testosteroni")
        & (df["comments"] == "10000 nM thiamine")
    ]
    df_oa = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
        & (df["comments"] == "10000 nM thiamine")
    ]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250328_ct_oa_thiamine_gradient/data/measurements.csv"
    )

    fig = go.Figure()
    for i, lg in enumerate(df_ct["linegroup"]):
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        slope = linregress(x[:36], np.log(y[:36]))[0]
        print("Slope Ct", slope)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="Ct",
                showlegend=False,
                line=dict(color=colors["ct"]),
            )
        )

    if model:
        p = parse_params()
        p["D"] = 0
        p["N01"] = y[0]
        p["q1_1"] = 0.028
        xs = data[lg + "_time"].to_numpy()
        Y_ct = odeint(ct_mono, [p["N01"], p["M1"]], xs, args=(p,))
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=Y_ct[:, 0],
                name="<i>C. testosteroni</i><br>model",
                marker=dict(color=colors["ct"]),
                mode="lines",
                line=dict(dash="dash"),
            )
        )

    for i, lg in enumerate(df_oa["linegroup"]):
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        slope = linregress(x[:36], np.log(y[:36]))[0]
        print("Slope Oa", slope)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="Oa 10 μM thiamine",
                showlegend=False,
                line=dict(color=colors["oa"]),
                mode="lines",
            )
        )

    if model:
        p["N02"] = y[0]
        Y_oa = odeint(oa_mono, [p["N02"], p["M1"]], xs, args=(p,))
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=Y_oa[:, 0],
                name="<i>O. anthropi</i><br>model",
                marker=dict(color=colors["oa"]),
                line=dict(dash="dash"),
                mode="lines",
            )
        )
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/metadata.csv"
    )
    df = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_oa_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
    ]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/measurements.csv"
    )
    for i, lg in enumerate(df[df["comments"] == "0 nM thiamine"]["linegroup"]):
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="Oa no thiamine",
                showlegend=False,
                line=dict(color=colors["oa"], dash="3px"),
                mode="lines",
            )
        )
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250328_oa_in_ct_OD_gradient/data/metadata.csv"
    )
    df = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_oa_in_spent_media_of_ct")
        & (df["species"] == "Ochrobactrum anthropi")
    ]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250328_oa_in_ct_OD_gradient/data/measurements.csv"
    )
    for i, lg in enumerate(df[df["comments"] == "0.37 OD of Ct"]["linegroup"]):
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="Oa in spent<br>media of Ct",
                showlegend=False,
                line=dict(color=colors["oa"], dash="1px"),
                mode="lines",
                marker=dict(color=colors["oa"]),
            )
        )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)],
            showgrid=True,
            zeroline=True,
            dtick=10,
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            # range=[0, 0.35],
            showgrid=True,
            # dtick=0.1,
            title="OD",
            ticks="inside",
        ),
        width=width,
        height=height,
        title="Batch cultures in minimal media",
    )
    fig = style_plot(
        fig,
        font_size=11,
        buttom_margin=20,
        top_margin=20,
        left_margin=20,
        right_margin=20,
    )
    fig.write_image("plots/experiments/growth_curves_ct_oa.svg")


def chemostat_ct_oa_community_cross_feeding():
    df = get_cfus()
    reactors = ["M0", "M1", "M2"]
    species = ["ct", "oa"]
    ct_oa = df[df["experiment"] == "ct_oa"]
    legend = {
        "ct": "Ct",
        "oa": "Oa",
        "ct_oa_thiamine": "A + T",
        "ct_oa": "A",
    }

    fig = go.Figure()
    for s in species:
        for i, r in enumerate(reactors):
            data = ct_oa[(ct_oa["reactor"] == r) & (ct_oa["species"] == s)]
            fig.add_trace(
                go.Scatter(
                    x=data["sample_time"],
                    y=data["average"],
                    error_y=dict(
                        type="data", array=data["stdev"].to_list(), visible=True
                    ),
                    name=legend[s],
                    showlegend=False,
                    mode="lines+markers",
                    marker=dict(
                        color=colors[s],
                        # line=dict(color="black", width=1.2),
                    ),
                    line=dict(color=colors[s]),
                    opacity=0.8,
                ),
            )
    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            # range=[0, 52], dtick=12),
            ticks="inside",
        ),
        yaxis=dict(
            title="CFUs/mL",
            type="log",
            range=[6, 9],
            dtick=1,
            exponentformat="power",
            showexponent="all",
            ticks="inside",
        ),
        width=150,
        height=150,
        title="Thiamine-free",
    )
    fig = style_plot(
        fig,
        font_size=11,
        right_margin=0,
        left_margin=45,
        buttom_margin=30,
        top_margin=20,
        marker_size=7,
        line_thickness=1.5,
    )

    fig.write_image("plots/experiments/chemostat_ct_oa_cross_feeding.svg")


def chemostat_ct_oa_thiamine():
    df = get_cfus()
    reactors = ["M0", "M1", "M2"]
    species = ["ct", "oa"]
    ct_oa = df[df["experiment"] == "ct_oa"]
    legend = {
        "ct": "Ct",
        "oa": "Oa",
        "ct_oa_thiamine": "A + T",
        "ct_oa": "A",
    }
    ct_oa_thiamine = df[df["experiment"] == "ct_oa_thiamine"]
    fig = go.Figure()
    for s in species:
        for i, r in enumerate(reactors):
            data = ct_oa_thiamine[
                (ct_oa["reactor"] == r) & (ct_oa_thiamine["species"] == s)
            ]
            fig.add_trace(
                go.Scatter(
                    x=data["sample_time"],
                    y=data["average"],
                    error_y=dict(
                        type="data", array=data["stdev"].to_list(), visible=True
                    ),
                    name=legend[s],
                    showlegend=False,
                    line=dict(color=colors[s]),
                    opacity=0.8,
                ),
            )
    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="CFUs/mL",
            type="log",
            range=[6, 9],
            dtick=1,
            exponentformat="power",
            showexponent="all",
            ticks="inside",
        ),
        width=width,
        height=height,
        title="Thiamine-added",
    )
    fig = style_plot(
        fig,
        font_size=11,
        right_margin=0,
        left_margin=45,
        buttom_margin=30,
        top_margin=20,
        marker_size=7,
    )
    fig.write_image("plots/experiments/chemostat_ct_oa_thiamine.svg")


def chemostat_ct_oa_stats():
    df = get_cfus()
    ct_oa = df[df["experiment"] == "ct_oa"]
    ct_oa_thiamine = df[df["experiment"] == "ct_oa_thiamine"]
    ct_oa_thiamine = ct_oa_thiamine[ct_oa_thiamine["sample_time"] != 0]
    ss = ct_oa_thiamine.loc[ct_oa_thiamine["sample_time"] != 0].copy()

    wide = ss.pivot_table(
        index=["reactor", "sample_time"],
        columns="species",
        values="average",
        aggfunc="mean",
    ).dropna(subset=["ct", "oa"])

    # log10(Oa/Ct) = log10(Oa) - log10(Ct)
    wide["d"] = np.log10(wide["oa"]) - np.log10(wide["ct"])
    dlong = wide.reset_index()  # columns: reactor, sample_time, ct, oa, d

    m = smf.mixedlm("d ~ 1", dlong, groups=dlong["reactor"])
    r = m.fit(reml=False)

    est = r.params["Intercept"]  # mean log10(Oa/Ct)
    se = r.bse["Intercept"]

    # one-sided p for Oa > Ct
    z = est / se
    p_one = 1 - stats.norm.cdf(z)

    # fold-change Oa/Ct with ~95% CI
    ci = (est - 1.96 * se, est + 1.96 * se)
    print("mean log10(Oa/Ct):", est, "one-sided p:", p_one)
    print("Oa/Ct fold:", 10**est, "CI:", (10 ** ci[0], 10 ** ci[1]))


def statistics_ct_oa_chemostat(eps=1e-12):
    df = get_cfus()
    ct_oa_thiamine = df[df["experiment"] == "ct_oa_thiamine"]
    # 1) filter time 0
    ss = ct_oa_thiamine.loc[ct_oa_thiamine["sample_time"] != 0].copy()

    # 2) wide table
    wide = ss.pivot_table(
        index=["reactor", "sample_time"],
        columns="species",
        values="average",
        aggfunc="mean",
    ).dropna(subset=["ct", "oa"])

    # 3) guard against non-positive values (log not defined)
    # if you *expect* zeros from detection limits, using eps is common:
    wide["ct_pos"] = wide["ct"].astype(float).clip(lower=eps)
    wide["oa_pos"] = wide["oa"].astype(float).clip(lower=eps)

    # 4) log10 ratio
    wide["d"] = np.log10(wide["oa_pos"]) - np.log10(wide["ct_pos"])
    dlong = wide.reset_index()

    # 5) mixed model: random intercept per reactor
    m = smf.mixedlm("d ~ 1", dlong, groups=dlong["reactor"])
    r = m.fit(reml=False)

    est = float(r.params["Intercept"])  # mean log10(Oa/Ct)
    se = float(r.bse["Intercept"])

    # Wald z-test for H0: est=0
    z = est / se

    # one-sided p for Oa > Ct (H1: est > 0)
    p_one = stats.norm.sf(z)  # same as 1 - cdf(z)

    # (optional) two-sided p
    p_two = 2 * stats.norm.sf(abs(z))

    # 95% CI on log10 scale, then back-transform to fold-change
    ci_log = (est - 1.96 * se, est + 1.96 * se)
    fold = 10**est
    fold_ci = (10 ** ci_log[0], 10 ** ci_log[1])

    out = {
        "mean_log10_Oa_over_Ct": est,
        "se": se,
        "z": z,
        "p_one_sided_Oa_gt_Ct": p_one,
        "p_two_sided": p_two,
        "fold_Oa_over_Ct": fold,
        "fold_CI95": fold_ci,
        "result": r,  # keep the fitted model if you want r.summary()
    }

    print(
        f"mean log10(Oa/Ct) = {out['mean_log10_Oa_over_Ct']:.3f} (SE {out['se']:.3f}), z={out['z']:.2f}"
    )
    print(
        f"one-sided p (Oa > Ct) = {out['p_one_sided_Oa_gt_Ct']:.3g}  |  two-sided p = {out['p_two_sided']:.3g}"
    )
    print(
        f"Oa/Ct fold-change = {out['fold_Oa_over_Ct']:.2f}  (95% CI {out['fold_CI95'][0]:.2f}–{out['fold_CI95'][1]:.2f})"
    )


def calibration_curve(cal: pd.DataFrame, title="OD calibration: OD = k·log10(R) + b"):
    """
    cal: DataFrame indexed by reactor (or with 'reactor' column),
         columns: OD0, R0, OD1, R1
    """
    cal = cal.copy()
    if "reactor" in cal.columns:
        cal = cal.set_index("reactor")

    # compute k and b per reactor (no sympy needed)
    cal["k"] = (cal["OD1"] - cal["OD0"]) / (np.log10(cal["R1"]) - np.log10(cal["R0"]))
    cal["b"] = cal["OD0"] - cal["k"] * np.log10(cal["R0"])

    fig = go.Figure()

    # choose an R-range for the line (across both points)
    R_all = np.r_[cal["R0"].to_numpy(float), cal["R1"].to_numpy(float)]
    R_min, R_max = float(np.nanmin(R_all)), float(np.nanmax(R_all))
    R_line = np.logspace(np.log10(R_min), np.log10(R_max), 200)
    names = ["Replicate 1", " Replicate 2", "Replicate 3"]
    for i, (reactor, row) in enumerate(cal.iterrows()):
        color = list(colors_metabolites.values())[i]

        # points
        fig.add_trace(
            go.Scatter(
                x=[row["R0"], row["R1"]],
                y=[row["OD0"], row["OD1"]],
                mode="markers",
                name=f"{names[i]}",
                marker=dict(size=9, color=color, line=dict(width=1)),
                showlegend=True,
            )
        )

        # line
        OD_line = row["k"] * np.log10(R_line) + row["b"]
        fig.add_trace(
            go.Scatter(
                x=R_line,
                y=OD_line,
                mode="lines",
                name=f"{names[i]} fit",
                line=dict(color=color, width=2),
                showlegend=False,
            )
        )

    fig.update_layout(
        title=title,
        xaxis=dict(title="Raw reading R", type="log", ticks="inside"),
        yaxis=dict(title="Calibrated OD", ticks="inside"),
        legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99),
        template="simple_white",
        width=270,
        height=180,
    )
    return fig


def plot_calibration_curve_ct_oa_mono():
    df = pd.read_csv("~/ChiBioFlow/data/at_oa/250320_ct_mono/calibration.csv")
    fig = calibration_curve(df)
    fig = style_plot(fig, font_size=11, marker_size=7, line_thickness=1.5)
    fig.write_image("plots/experiments/od_calibration_curve_ct_mono.svg")

    df = pd.read_csv("~/ChiBioFlow/data/at_oa/250310_oa_mono/calibration.csv")
    fig = calibration_curve(df)
    fig = style_plot(fig, font_size=11, marker_size=7, line_thickness=1.5)
    fig.write_image("plots/experiments/od_calibration_curve_oa_mono.svg")


def plot_ct_in_spent_chemostat_media():
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/metadata.csv"
    )
    df_ct = df[df["species"] == "Comamonas testosteroni"]
    df_ct = df_ct.sort_values(by="exp_ID")
    df_oa = df[df["species"] == "Ochrobactrum anthropi"]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/measurements.csv"
    )
    figs = [go.Figure(), go.Figure()]
    ys = []
    sources = [
        "Spent media Ct",
        "Spent media Oa",
    ]
    cvs = []
    for j, cs in enumerate(sources):
        line_counter = {"Batch 1": 0, "Batch 2": 0, "Batch 3": 0}
        for i, (lg, comment) in enumerate(
            df_ct[df_ct["carbon_source"] == cs][["linegroup", "comments"]].values
        ):
            if line_counter[comment] == 0:
                cv = pd.DataFrame(columns=["Time", "OD", "Well", "Strain"])
                x = data[lg + "_time"]
                y = data[lg + "_measurement"]
                ys.append(y.to_numpy())
                figs[j].add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        name=cs,
                        legendgroup=cs,
                        line=dict(color=colors["ct"]),
                        showlegend=(i == 0),
                        opacity=(1 if j == 0 else 1),
                    )
                )
                line_counter[comment] += 1
                cv["Time"], cv["OD"], cv["Well"], cv["Strain"] = x, y, lg, cs
                cv = cv[cv["Time"] > 18]
                cvs.append(cv)

    cv = pd.concat(cvs)
    cv_ct_oa = cv[cv["Strain"] == "Spent media Oa"]

    m = curveball.models.fit_model(
        cv_ct_oa,
        PLOT=False,
        PRINT=False,
        param_guess={"y0": 0.009},
        param_fix=["y0"],
    )
    m = sorted(m, key=lambda r: r.model.name)
    figs[1].add_trace(
        go.Scatter(
            x=m[1].userkws["t"],
            y=m[1].best_fit,
            mode="lines",
            line=dict(dash="1px 1px", color="black"),
            name="fit",
            showlegend=False,
        )
    )

    figs[1].update_layout(
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD", ticks="inside"),
        showlegend=False,
        width=175,
        height=160,
        title=f"Ct in spent chemo-<br>stat media of Oa r: {m[1].params['r'].value:.2f}",
    )
    figs[1] = style_plot(
        figs[1],
        font_size=11,
        buttom_margin=25,
        left_margin=35,
        right_margin=10,
        top_margin=30,
    )
    print(f"Ct in spent chemo-stat media of Oa r: {m[1].params['r'].value:.2f}")
    figs[1].write_image("plots/experiments/ct_grown_in_oa_chemostat.svg")

    cv_ct_ct = cv[cv["Strain"] == "Spent media Ct"]

    m = curveball.models.fit_model(
        cv_ct_ct,
        PLOT=False,
        PRINT=False,
        param_guess={"y0": 0.009},
        param_fix=["y0"],
    )
    m = sorted(m, key=lambda r: r.model.name)
    figs[0].add_trace(
        go.Scatter(
            x=m[1].userkws["t"],
            y=m[1].best_fit,
            mode="lines",
            line=dict(dash="1px 1px", color="black"),
            name="fit",
            showlegend=False,
        )
    )

    figs[0].update_layout(
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD", ticks="inside"),
        showlegend=False,
        width=175,
        height=160,
        title=f"Ct in spent chemo-<br>stat media of Ct r: {m[1].params['r'].value:.2f}",
    )
    figs[0] = style_plot(
        figs[0],
        font_size=11,
        buttom_margin=25,
        left_margin=35,
        right_margin=10,
        top_margin=30,
    )
    figs[0].write_image("plots/experiments/ct_grown_in_ct_chemostat.svg")
    print(f"Ct in spent chemo-stat media of Ct r: {m[1].params['r'].value:.2f}")


def plot_oa_in_spent_chemostat_media():
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/metadata.csv"
    )
    df_ct = df[df["species"] == "Comamonas testosteroni"]
    df_ct = df_ct.sort_values(by="exp_ID")
    df_oa = df[df["species"] == "Ochrobactrum anthropi"]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/measurements.csv"
    )
    figs = [go.Figure(), go.Figure()]
    css = ["Spent media Ct", "Spent media Oa"]
    for j, cs in enumerate(css):
        line_counter = {"Batch 1": 0, "Batch 2": 0, "Batch 3": 0}
        for i, (lg, comment) in enumerate(
            df_oa[df_oa["carbon_source"] == cs][["linegroup", "comments"]].values
        ):
            if line_counter[comment] == 0:
                x = data[lg + "_time"]
                y = data[lg + "_measurement"]
                figs[j].add_trace(
                    go.Scatter(
                        x=x[5:],
                        y=y[5:],
                        name=cs,
                        legendgroup=cs,
                        mode="lines",
                        line=dict(
                            color=colors["oa"],
                        ),
                        showlegend=False,
                    )
                )
                line_counter[comment] += 1
    titles = [
        "Oa in spent chemo-<br>stat media of Ct",
        "Oa in spent chemo-<br>stat media of Oa",
    ]
    for i, title in enumerate(titles):
        figs[i].update_layout(
            xaxis=dict(title="Time [h]", range=[0, 72], dtick=12),
            yaxis=dict(title="OD", range=[0, 0.03], dtick=0.01),
            title=title,
            width=175,
            height=160,
            showlegend=False,
        )
        figs[i] = style_plot(
            figs[i],
            font_size=11,
            buttom_margin=25,
            left_margin=35,
            right_margin=10,
            top_margin=30,
        )
        figs[i].write_image(
            "plots/experiments/oa_grown_in_{}_chemostat.svg".format(
                css[i].lower().replace(" ", "_")
            )
        )


def oa_mono_no_thiamine():
    Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/250310_oa_mono")
    df = calibration_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250310_oa_mono/calibration.csv", Ms
    )
    df = df[df["reactor"] == "M3"]

    fig = go.Figure()

    # experimental data
    fig.add_trace(
        go.Scatter(
            x=df["exp_time"][4:-10],
            y=df["od_calibrated"][4:-10],
            mode="lines+markers",
            name="Chemostat",
            showlegend=False,
            marker=dict(color=colors["oa"]),
            line=dict(color=colors["oa"]),
        )
    )

    # washout reference: starts at OD 0.15 after 20 h
    D = 0.15
    t0 = 20
    od0 = 0.15

    t_wash = np.linspace(t0, df["exp_time"].max(), 300)
    od_wash = od0 * np.exp(-D * (t_wash - t0))

    fig.add_trace(
        go.Scatter(
            x=t_wash,
            y=od_wash,
            mode="lines",
            name="Washout",
            showlegend=False,
            line=dict(
                color=colors["oa"],
                dash="dash",
                width=2,
            ),
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(title="OD600", ticks="inside"),
        title="Oa in chemostat without thiamine",
        width=width,
        height=height,
        showlegend=False,
    )

    fig = style_plot(fig, font_size=11, left_margin=20, buttom_margin=20)
    fig.write_image("plots/experiments/oa_mono_no_thiamine.svg")


oa_mono_no_thiamine()


def oa_thiamine_gradient():
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
    keep = list(colors.keys())[:6]

    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/metadata.csv"
    )
    df = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_oa_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
    ]
    mask = []
    keep = list(colors.keys())[:6]
    for c in df["comments"]:
        if c in keep:
            mask.append(True)
        else:
            mask.append(False)
    df = df[mask]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/measurements.csv"
    )
    fig = go.Figure()

    for i, (lg, t_conc) in enumerate(zip(df["linegroup"], df["comments"])):
        x = data[lg + "_time"]
        y = data[lg + "_measurement"]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name=" ".join(t_conc.split(" ")[:2]),
                showlegend=(i <= 5),
                marker=dict(color=colors[t_conc]),
            )
        )

    p = parse_params()
    p["D"] = 0
    xs = np.linspace(0, max(x), 200)
    Y = odeint(
        thiamine_supply,
        [0, 0.01, 7.5, 1],
        xs,
        args=(p,),
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y[:, 1],
            name="Model",
            line=dict(dash="dot", color="black"),
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)],
            showgrid=True,
            zeroline=True,
            dtick=12,
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            range=[0, 0.4],
            showgrid=True,
            zeroline=True,
            dtick=0.1,
            title="OD",
            ticks="inside",
        ),
        width=150,
        height=150,
        showlegend=False,
        title="Oa thiamine gradient",
    )
    fig = style_plot(
        fig,
        font_size=11,
        left_margin=30,
        right_margin=0,
        buttom_margin=30,
        top_margin=20,
    )
    fig.write_image("plots/experiments/oa_thiamine_gradient.svg")


def ct_oa_no_cs():
    legend = {
        "ct": "Ct",
        "oa": "Oa",
        "ct_oa_thiamine": "A + T",
        "ct_oa": "A",
    }
    cfus = get_cfus()
    cfus = cfus[(cfus["experiment"] == "no_cs") & (cfus["species"].isin(["ct", "oa"]))]
    fig = go.Figure()
    for i, s in enumerate(cfus["species"].unique()):
        df = cfus[cfus["species"] == s]
        fig.add_trace(
            go.Scatter(
                x=df["sample_time"],
                y=df["average"],
                error_y=dict(type="data", array=df["stdev"].to_list(), visible=True),
                name=legend[s],
                showlegend=False,
                line=dict(color=colors[s]),
            ),
        )
    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            # range=[0, 42],
            # dtick=12),
            ticks="inside",
        ),
        yaxis=dict(
            title="CFUs/mL",
            type="log",
            range=[5, 9],
            ticks="inside",
        ),
        width=150,
        height=180,
        title="No carbon source",
    )
    fig = style_plot(
        fig,
        font_size=11,
        right_margin=0,
        left_margin=45,
        buttom_margin=30,
        top_margin=20,
        marker_size=7,
    )
    fig.write_image("plots/experiments/no_cs.svg")
