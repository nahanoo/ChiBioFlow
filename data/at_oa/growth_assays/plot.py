import pandas as pd
import plotly.express as px
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
import curveball
from sympy import symbols, solve, diff, Matrix, Eq
from sympy import init_printing
from os.path import join

init_printing()


def dump_dfs():
    at = pd.read_csv("at.tsv", sep="\t")
    at.insert(len(at.columns), "name", "at")
    oa = pd.read_csv("oa.tsv", sep="\t")
    oa.insert(len(oa.columns), "name", "oa")

    df = pd.concat([at, oa])
    dfs = []

    for i, c in enumerate(df.columns[1:-1]):
        if i % 3 == 0:
            carbon_source = i
        tmp = df[["Time", c, "name"]]
        out = pd.DataFrame()
        out["time"], out["OD"], out["well"], out["strain"], out["carbon_source"] = (
            tmp["Time"],
            tmp[c],
            c,
            tmp["name"],
            carbon_source,
        )
        dfs.append(out)

    df = pd.concat(dfs)
    df.to_csv("carbon_sources.csv", index=False)
    fum = df[df["carbon_source"] == 21]
    ala = df[df["carbon_source"] == 36]
    glut = df[df["carbon_source" == 63]]

    at = fum[fum["strain"] == "at"]
    at_rates = []
    for well in set(at["well"]):
        OD = at[at["well"] == well]["OD"]
        at_rates += list(np.gradient(OD) / OD * 2)
    at.insert(len(at.columns), "growth_rates", at_rates)

    oa = fum[fum["strain"] == "oa"]
    oa_rates = []
    for well in set(oa["well"]):
        OD = oa[oa["well"] == well]["OD"]
        oa_rates += list(np.gradient(OD) / OD * 2)
    oa.insert(len(oa.columns), "growth_rates", oa_rates)
    pd.concat([at, oa]).to_csv("fumarate.csv", index=False)

    at = ala[ala["strain"] == "at"]
    at_rates = []
    for well in set(at["well"]):
        OD = at[at["well"] == well]["OD"]
        at_rates += list(np.gradient(OD) / OD * 2)
    at.insert(len(at.columns), "growth_rates", at_rates)

    oa = ala[ala["strain"] == "oa"]
    oa_rates = []
    for well in set(oa["well"]):
        OD = oa[oa["well"] == well]["OD"]
        oa_rates += list(np.gradient(OD) / OD * 2)
    oa.insert(len(oa.columns), "growth_rates", oa_rates)
    pd.concat([at, oa]).to_csv("alanine.csv", index=False)

    at = glut[glut["strain"] == "at"]
    at_rates = []
    for well in set(at["well"]):
        OD = at[at["well"] == well]["OD"]
        at_rates += list(np.gradient(OD) / OD * 2)
    at.insert(len(at.columns), "growth_rates", at_rates)

    oa = glut[glut["strain"] == "oa"]
    oa_rates = []
    for well in set(oa["well"]):
        OD = oa[oa["well"] == well]["OD"]
        oa_rates += list(np.gradient(OD) / OD * 2)
    oa.insert(len(oa.columns), "growth_rates", oa_rates)
    pd.concat([at, oa]).to_csv("alanine.csv", index=False)


dir = join("/", "home", "eric", "ChiBioFlow", "data", "at_oa", "growth_assays")
fum = pd.read_csv(join(dir, "fumarate.csv"))
ala = pd.read_csv(join(dir, "alanine.csv"))
css = pd.read_csv(join(dir, "carbon_sources.csv"))


def plot_carbon_sources():
    return px.line(
        css,
        x="time",
        y="OD",
        facet_col="carbon_source",
        color="strain",
        facet_col_wrap=6,
        line_group="well",
    )


def plot_fumarate():
    px.line(fum, x="time", y="OD", color="strain", line_group="well").show()
    px.line(fum, x="time", y="growth_rates", color="strain", line_group="well").show()


def plot_alanine():
    OD = px.line(ala, x="time", y="OD", color="strain", line_group="well")
    growth_rates = px.line(
        ala, x="time", y="growth_rates", color="strain", line_group="well"
    )
    return OD, growth_rates


def plot_malatat():
    mal = css[css["carbon_source"] == 63]
    px.line(mal, x="time", y="OD", color="strain", line_group="well").show()


def get_df(well, strain):
    df = css[(css["well"] == well) & (css["strain"] == strain)]
    df.columns = ["Time", "OD", "well", "strain", "carbon_source"]
    df = df.dropna()
    df["Time"] = np.linspace(0.5, 75.5, len(df))
    return df


def fit(well, strain):
    df = get_df(well, strain)
    models, fig, ax = curveball.models.fit_model(
        df, PLOT=True, PRINT=False, param_fix=["y0"]
    )
    return fig


def growth_paramters():
    N, Y, M, K, D, r, S = symbols("N,Y,M,K,D,r,S")

    fX = Eq(N, Y * (M - K * D / (r - D)))
    fK = solve(fX, K)[0]
    K_value = fK.subs({"N": 0.13, "Y": 0.03, "M": 30, "r": 0.34, "D": 0.16})

    fY = Eq(N / (M - S), Y)
    fS = solve(fY, S)[0]
    ss = fS.subs({"N": 0.19, "Y": 0.9 / 30, "M": 30})
    print(K_value)


def regression_oa():
    time = np.array([0.0, 0.5, 1.0, 1.5, 2.0])[:2]
    od = np.array([0.188, 0.248, 0.316, 0.38, 0.429])[:2]

    # Log transform the OD data
    log_od = np.log(od)

    # Fit a linear regression to the data points
    slope, intercept, r_value, p_value, std_err = linregress(time, log_od)

    # Plotting to visualize
    plt.plot(time, log_od, "o", label="Log(OD)")
    plt.plot(time, intercept + slope * time, "r", label="Linear Fit")
    plt.xlabel("Time")
    plt.ylabel("Log(OD)")
    plt.title("Exponential Growth Phase")
    plt.legend()
    plt.show()

    print(f"The estimated growth rate (slope) is: {slope}")


def regression_at():
    # Example data (time, OD)

    # Extract time and OD columns
    od = np.array(at.loc[:60]["OD"])[21:23]
    time = np.linspace(0, len(od) * 0.5, len(od))

    # Log transform the OD data
    log_od = np.log(od)

    # Fit a linear regression to the data points
    slope, intercept, r_value, p_value, std_err = linregress(time, log_od)

    # Plotting to visualize
    plt.plot(time, log_od, "o", label="Log(OD)")
    plt.plot(time, intercept + slope * time, "r", label="Linear Fit")
    plt.xlabel("Time")
    plt.ylabel("Log(OD)")
    plt.title("Exponential Growth Phase")
    plt.legend()
    plt.show()

    print(f"The estimated growth rate (slope) is: {slope}")
