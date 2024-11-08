import pandas as pd
import numpy as np
import scipy.optimize as optimize
from os.path import join
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.stats import linregress
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from scipy.integrate import odeint
import scipy.optimize as optimize
from scipy.signal import savgol_filter

species_dict = {"Comamonas testosteroni": "ct", "Ochrobactrum anthropi": "oa"}
colors = {
    "Comamonas testosteroni": "#7570B3",
    "Ochrobactrum anthropi": "#D95F02",
    "total": "#565e91",
}


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
        [0.0],
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


def parse_data():
    df = pd.DataFrame(
        columns=["species", "concentration", "time", "OD", "v", "q", "Km"])
    meta_full = pd.read_csv(
        "/home/eric/curves/metadata/pooled_df_joint_metadata.csv")
    raw = pd.read_csv(
        "/home/eric/curves/export/240921_carbon_source_gradients/measurement_data.csv"
    )
    raw.index = raw["linegroup"]
    for species in ["Comamonas testosteroni", "Ochrobactrum anthropi"]:
        mask = (meta_full["project"] == "240921_carbon_source_gradients") & (
            meta_full["species"] == species)
        meta = meta_full[mask]

        row, col = 0, 0
        fig = make_subplots(
            rows=3,
            cols=3,
            subplot_titles=list(set(meta["cs_conc"])),
            shared_yaxes=True,
            shared_xaxes=True,
        )
        for i, c in enumerate(set(meta["cs_conc"])):
            xs, ys = [], []
            c_sub = meta[meta["cs_conc"] == c]
            if i % 3 == 0:
                col = 1
                row += 1
            else:
                col += 1
            for l in c_sub["linegroup"]:
                x, y = (
                    raw.loc[l, "time"].to_numpy(),
                    raw.loc[l, "measurement"].to_numpy(),
                )

                ys.append(y), xs.append(x)
            df.loc[species.replace(" ", "_") + "_" + str(c)] = [
                species,
                c,
                np.average(np.array(xs), axis=0),
                np.average(np.array(ys), axis=0),
                None,
                None,
                None,
            ]
    return df


def fit_params(df):
    for species in set(df["species"]):
        for c in set(df[df["species"] == species]["concentration"]):
            for i, row in df[(df["species"] == species)
                             & (df["concentration"] == c)].iterrows():
                x, y = row["time"], row["OD"]
                if species == "Comamonas testosteroni":
                    q = 0.4 / 7.5
                    exp_phase = 6
                    stat_phase = 12
                if species == "Ochrobactrum anthropi":
                    exp_phase = 12
                    q = 0.5 / 7.5
                    stat_phase = 40
                x_window = x[x < exp_phase]
                y_window = y[:len(x_window)]
                v = linregress(x_window, np.log(y_window))[0]
                df.loc[i, "q"] = q
                df.loc[i, "v"] = v
                x_stat = x[x < stat_phase]
                y_stat = y[:len(x_stat)]
                Km = get_Km(x_stat, y_stat, c, 0.005, v, df.loc[i, "q"])
                df.loc[i, "Km"] = Km
    return df


def plot_plate_reader(fit=False):
    df = parse_data()
    figs = []
    if fit:
        df = fit_params(df)
    df = df.sort_values("concentration")
    for species in ["Comamonas testosteroni", "Ochrobactrum anthropi"]:
        row_f, col_f = 0, 0
        fig = make_subplots(
            rows=3,
            cols=3,
            subplot_titles=[
                str(c) + " mM Acetate"
                for c in sorted(list(set(df["concentration"])))
            ],
            horizontal_spacing=0.03,
            vertical_spacing=0.05,
            shared_yaxes=True,
            shared_xaxes=True,
        )
        if fit:
            showlegend = True
        else:
            showlegend = False
        for j, (i, row) in enumerate(df[df["species"] == species].iterrows()):
            if j % 3 == 0:
                col_f = 1
                row_f += 1
            else:
                col_f += 1
            fig.add_trace(
                go.Scatter(
                    x=row["time"],
                    y=row["OD"],
                    name="OD",
                    showlegend=showlegend,
                    line=dict(color=colors[species]),
                ),
                col=col_f,
                row=row_f,
            )
            if fit:
                sim = odeint(
                    monod,
                    [0.005, row["concentration"]],
                    row["time"],
                    args=(row["v"], row["Km"], row["q"]),
                )[:, 0]
                fig.add_trace(
                    go.Scatter(
                        x=row["time"],
                        y=sim,
                        line=dict(color=colors[species], dash="dash"),
                        name="Fit",
                        showlegend=showlegend,
                    ),
                    col=col_f,
                    row=row_f,
                )
            showlegend = False
        fig.update_layout(title=species)
        fig.update_xaxes(title="Time [h]", row=3, col=2)
        fig.update_yaxes(title="OD", row=2, col=1)
        figs.append(fig)
    return figs


def plate_reader_growth_rates():
    figs = []
    df = parse_data()
    df = fit_params(df)
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=["Comamonas testosteroni", "Ochrobactrum anthropi"],
        shared_yaxes=True,
    )
    for i, species in enumerate(
        ["Comamonas testosteroni", "Ochrobactrum anthropi"]):
        df_sub = df[df["species"] == species]
        df_sub = df_sub.sort_values("concentration")
        fig.add_trace(
            go.Scatter(
                x=df_sub["concentration"],
                y=df_sub["v"],
                line=dict(color=colors[species]),
                showlegend=False,
            ),
            col=i + 1,
            row=1,
        )
        fig.update_xaxes(title="Concentration [mM]", col=1)
        fig.update_yaxes(title="Max. growth rate [1/h]")
    figs.append(fig)
    return figs


def plate_reader_max_OD():
    figs = []
    meta_full = pd.read_csv(
        "/home/eric/curves/metadata/pooled_df_joint_metadata.csv")
    raw = pd.read_csv(
        "/home/eric/curves/export/240921_carbon_source_gradients/measurement_data.csv"
    )
    raw.index = raw["linegroup"]

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=["Comamonas testosteroni", "Ochrobactrum anthropi"],
        shared_yaxes=True,
    )
    yield_od = pd.DataFrame(
        columns=["species", "concentration", "avg_OD", "std_OD"])
    for j, species in enumerate(
        ["Comamonas testosteroni", "Ochrobactrum anthropi"]):
        mask = (meta_full["project"] == "240921_carbon_source_gradients") & (
            meta_full["species"] == species)
        meta = meta_full[mask]
        for i, c in enumerate(set(meta["cs_conc"])):
            ys = []
            c_sub = meta[meta["cs_conc"] == c]
            for l in c_sub["linegroup"]:
                x, y = (
                    raw.loc[l, "time"].to_numpy(),
                    raw.loc[l, "measurement"].to_numpy(),
                )
                ys.append(max(y))
            yield_od.loc[len(yield_od)] = [
                species, c, np.average(ys),
                np.std(ys)
            ]
        yield_od_sub = yield_od[yield_od["species"] == species]
        yield_od_sub = yield_od_sub.sort_values("concentration")
        fig.add_trace(
            go.Scatter(
                x=yield_od_sub["concentration"],
                y=yield_od_sub["avg_OD"],
                line=dict(color=colors[species]),
                showlegend=False,
            ),
            col=j + 1,
            row=1,
        )
        fig.update_xaxes(title="Concentration [mM]", col=1)
        fig.update_yaxes(title="Max. OD")
    figs.append(fig)
    return figs


def max_growth_rate_oa():
    df = parse_data()
    oa = df.loc['Ochrobactrum_anthropi_7.5']
    x, y = oa['time'], oa['OD']
    x = x[x < 7.5]
    y = y[:len(x)]
    v, b = linregress(x, np.log(y))[:2]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=np.log(y)))
    fig.add_trace(
        go.Scatter(x=x, y=[v * t + b for t in x], line=dict(dash='dash')))
    print(v)
    fig.show()


def max_growth_rate_ct():
    pass


df = parse_data()
ct = df.loc['Comamonas_testosteroni_7.5']
x, y = ct['time'], ct['OD']
x = x[x < 8.5]
y = y[:len(x)]
x, y = x[24:], y[24:]
v, b = linregress(x, np.log(y))[:2]
fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=np.log(y)))
fig.add_trace(go.Scatter(x=x, y=[v * t + b for t in x],
                         line=dict(dash='dash')))
print(v)
fig.show()
