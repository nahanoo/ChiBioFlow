import pandas as pd
from data_parser import get_cfus, get_od_chemostats, ct_oa_plate_reader
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from chibio_parser import *
from models import *
from scipy.stats import linregress
from style import *


def plot_chemostat_community():
    df = get_cfus()
    reactors = ["M0", "M1", "M2"]
    species = ["ct", "oa"]
    ct_oa = df[df["experiment"] == "ct_oa"]
    ct_oa_thiamine = df[df["experiment"] == "ct_oa_thiamine"]
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
                    line=dict(color=colors[s]),
                ),
            )
    fig.update_layout(
        xaxis=dict(title="Time [h]", range=[0, 52], dtick=12),
        yaxis=dict(
            title="CFUs/mL",
            type="log",
            range=[7, 9],
        ),
        width=150,
        height=150,
    )
    fig = style_plot(
        fig,
        font_size=11,
        line_thickness=1,
        right_margin=0,
        left_margin=45,
        buttom_margin=30,
        top_margin=0,
    )
    fig.write_image("plots/experiments/ct_oa_a.svg")

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
                ),
            )
    fig.update_layout(
        xaxis=dict(title="Time [h]", range=[0, 42], dtick=12),
        yaxis=dict(
            title="CFUs/mL",
            type="log",
            range=[7, 9],
        ),
        width=width,
        height=height,
    )
    fig = style_plot(
        fig,
        font_size=11,
        line_thickness=1,
        right_margin=0,
        left_margin=45,
        buttom_margin=30,
        top_margin=0,
    )
    fig.write_image("plots/experiments/ct_oa_a_t.svg")


def plot_oa_mono_chemostats():
    fig = go.Figure()

    df = get_od_chemostats()
    df = df[df["experiment"] == "oa_mono"]
    Ms = [df[df["reactor"] == M] for M in ["M0", "M1"]]
    df = get_od_chemostats()
    df = df[df["experiment"] == "oa_mono_repeat"]
    df["reactor"] = "M2"
    Ms.insert(2, df)
    for i, M in enumerate(Ms):
        x = M["exp_time"].to_numpy()
        y = M["od_calibrated"].to_numpy()
        fig.add_trace(
            go.Scatter(
                x=x[4:-10],
                y=y[4:-10],
                # name=M.loc[0, "reactor"],
                name="Chemostat",
                showlegend=True,
                marker=dict(color=colors["oa"]),
            )
        )

    fig.update_layout(
        xaxis=dict(
            range=[0, max(M["exp_time"])],
            showgrid=True,
            zeroline=True,
            dtick=12,
            title="Time [h]",
        ),
        yaxis=dict(range=[0, 0.5], dtick=0.1),
        title="<i>O. anthropi</i>",
        width=width,
        height=height,
        showlegend=False,
    )
    p = parse_params()

    p["N02"] = 0.1
    p["q2_1"] = 0.053
    Y = odeint(oa_mono, [p["N02"], p["M1"]], M["exp_time"], args=(p,))
    fig.add_trace(
        go.Scatter(
            x=M["exp_time"],
            y=Y[:, 0],
            name="Model",
            line=dict(color="black", dash="dot"),
            mode="lines",
        ),
    )
    fig = style_plot(
        fig, line_thickness=1, font_size=11, left_margin=20, buttom_margin=20
    )
    fig.write_image("plots/experiments/oa_mono.svg")
    df = get_od_chemostats()
    df = df[df["experiment"] == "oa_mono"]
    Ms = [df[df["reactor"] == M] for M in ["M3"]]
    fig = go.Figure()
    for i, M in enumerate(Ms):
        x = M["exp_time"].to_numpy()
        y = M["od_calibrated"].to_numpy()
        fig.add_trace(
            go.Scatter(
                x=x[4:-10],
                y=y[4:-10],
                # name=M.loc[0, "reactor"],
                name="Chemostat",
                showlegend=True,
                marker=dict(color=colors["oa"]),
            )
        )

    fig.update_layout(
        xaxis=dict(
            range=[0, max(M["exp_time"])],
            showgrid=True,
            zeroline=True,
            dtick=12,
            title="Time [h]",
        ),
        yaxis=dict(range=[0, 0.5], dtick=0.1),
        title="<i>O. anthropi</i>",
        width=width,
        height=height,
        showlegend=False,
    )
    fig = style_plot(
        fig, line_thickness=1, font_size=11, left_margin=20, buttom_margin=20
    )
    fig.write_image("plots/experiments/oa_mono_no_thiamine.svg")


def plot_oa_ct_plate_reader():
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
                marker=dict(color=colors["ct"]),
            )
        )
        p = parse_params()
    p["D"] = 0
    p["N01"] = y[0]
    p["q1_1"] = 0.028
    xs = data[lg + "_time"].to_numpy()
    Y_ct = odeint(ct_mono, [p["N01"], p["M1"]], xs, args=(p,))
    if model:
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
    """df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/metadata.csv"
    )
    df_oa = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_oa_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
        & (df["comments"] == "10000 nM thiamine")
    ]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/measurements.csv"
    )"""
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
    p["N02"] = y[0]
    Y_oa = odeint(oa_mono, [p["N02"], p["M1"]], xs, args=(p,))
    if model:
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
                line=dict(color=colors["oa"], dash="dash"),
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
                line=dict(color=colors["oa"], dash="dot"),
                mode="lines",
            )
        )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)], showgrid=True, zeroline=True, dtick=12, title="Time [h]"
        ),
        yaxis=dict(
            range=[0, 0.35],
            showgrid=True,
            dtick=0.1,
            title="OD",
        ),
        width=width,
        height=width,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        font_size=11,
        buttom_margin=10,
        top_margin=10,
        left_margin=10,
        right_margin=10,
    )
    fig.write_image("plots/experiments/ct_oa_plate_reader.svg")


def ct_oa_plate_reader_fit():
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
            line=dict(dash="dot"),
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
                line=dict(
                    color=colors["oa"],
                ),
                mode="lines",
            )
        )
    p["N02"] = y[0]
    Y_oa = odeint(oa_mono, [p["N02"], p["M1"]], xs, args=(p,))
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=Y_oa[:, 0],
            name="<i>O. anthropi</i><br>model",
            marker=dict(color=colors["oa"]),
            line=dict(dash="dot"),
            mode="lines",
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)], showgrid=True, zeroline=True, dtick=12, title="Time [h]"
        ),
        yaxis=dict(
            range=[0, 0.35],
            showgrid=True,
            dtick=0.1,
            title="OD",
        ),
        width=width,
        height=width,
        showlegend=False,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        font_size=11,
        buttom_margin=10,
        top_margin=10,
        left_margin=10,
        right_margin=10,
    )
    fig.write_image("plots/experiments/ct_oa_plate_reader_fit.svg")


def ct_oa_max_growth_rate():
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
    p = parse_params()
    p["D"] = 0
    xs = np.linspace(0, 72, 200)
    x = np.average([data[lg + "_time"] for lg in df_ct["linegroup"]], axis=0)
    y = np.average([data[lg + "_measurement"] for lg in df_ct["linegroup"]], axis=0)
    x = x[x < 4]
    y = y[: len(x)]
    slope, intercept, r_value, p_value, std_err = linregress(x, np.log(y))
    (print("Slope Ct", slope))
    fit = [slope * i + intercept for i in x]
    fig = go.Figure()
    for i, lg in enumerate(df_ct["linegroup"]):
        x = data[lg + "_time"]
        y = data[lg + "_measurement"]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=np.log(y),
                name="<i>C. testosteroni</i>",
                legendgroup="<i>C. testosteroni</i>",
                showlegend=(i == 0),
                marker=dict(color=colors["ct"]),
            )
        )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=fit,
            name="fit",
            line=dict(dash="dash", color="black"),
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)], showgrid=True, zeroline=True, dtick=6, title="Time [h]"
        ),
        yaxis=dict(
            range=[-4, 0], showgrid=True, zeroline=True, dtick=0.5, title="log(OD)"
        ),
    )
    x = np.average([data[lg + "_time"] for lg in df_oa["linegroup"]], axis=0)
    y = np.average([data[lg + "_measurement"] for lg in df_oa["linegroup"]], axis=0)
    x = x[x < 4]
    y = y[: len(x)]
    slope, intercept, r_value, p_value, std_err = linregress(x, np.log(y))
    print("Slope Oa", slope)
    fit = [p["v2_1"] * i + np.log(p["N02"]) for i in x]
    for i, lg in enumerate(df_oa["linegroup"]):
        x = data[lg + "_time"]
        y = data[lg + "_measurement"]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=np.log(y),
                name="<i>O. anthropi</i>",
                legendgroup="<i>O. anthropi</i>",
                showlegend=(i == 0),
                marker=dict(color=colors["oa"]),
            )
        )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=fit,
            name="fit",
            line=dict(dash="dash", color="black"),
            showlegend=False,
        )
    )
    fig = style_plot(fig, line_thickness=1.7)
    fig.write_image("plots/experiments/ct_max_growth_rate.svg")


def plot_ct_mono_chemostats():
    fig = go.Figure()
    df = get_od_chemostats()
    df = df[df["experiment"] == "ct_mono"]
    Ms = [df[df["reactor"] == M] for M in sorted(set(df["reactor"]))]
    for i, M in enumerate(Ms):
        x = M["exp_time"].to_numpy()
        y = M["od_calibrated"].to_numpy()
        fig.add_trace(
            go.Scatter(
                x=x[4:-10],
                y=y[4:-10],
                # name=M.loc[0, "reactor"],
                name="Chemostat",
                showlegend=(True if i == 0 else False),
                marker=dict(color="#7570B3"),
            )
        )

    fig.update_layout(
        xaxis=dict(
            range=[0, max(M["exp_time"])],
            showgrid=True,
            zeroline=True,
            dtick=6,
            title="Time [h]",
        ),
        yaxis=dict(range=[0, 0.5], dtick=0.1, title="OD"),
        title="<i>C. testosteroni</i>",
        width=width,
        height=height,
        showlegend=False,
    )
    p = parse_params()

    p["N01"] = 0.05
    p["q1_1"] = 0.047
    Y = odeint(ct_mono, [p["N01"], p["M1"]], M["exp_time"], args=(p,))
    fig.add_trace(
        go.Scatter(
            x=M["exp_time"],
            y=Y[:, 0],
            name="Model",
            line=dict(color="black", dash="dot"),
            mode="lines",
        ),
    )
    fig = style_plot(
        fig, line_thickness=1, font_size=11, left_margin=10, buttom_margin=20
    )
    fig.write_image("plots/experiments/ct_mono.svg")
    colors = ["#6A5ACD", "#7570B3"]
    fig = go.Figure()
    df = get_od_chemostats()
    df = df[df["experiment"] == "ct_mono_old"]
    df = df[df["exp_time"] <= 24]
    Ms = [df[df["reactor"] == M] for M in sorted(set(df["reactor"]))]
    for i, M in enumerate(Ms):
        x = M["exp_time"].to_numpy()
        y = M["od_calibrated"].to_numpy()
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name=M.loc[0, "reactor"],
                showlegend=True,
                marker=dict(color=colors[i]),
            )
        )

    fig.update_layout(
        xaxis=dict(
            range=[0, max(M["exp_time"])],
            showgrid=True,
            zeroline=True,
            dtick=6,
            title="Time [h]",
        ),
        yaxis=dict(
            range=[0, 0.5], showgrid=True, zeroline=True, dtick=0.05, title="log(OD)"
        ),
        title="<i>O. anthropi</i>",
        # width=width,
        # height=height,
    )
    p = parse_params()

    p["N02"] = 0.03
    p["q1_1"] = 0.022
    Y = odeint(ct_mono, [p["N02"], p["M1"]], M["exp_time"], args=(p,))
    fig.add_trace(
        go.Scatter(
            x=M["exp_time"],
            y=Y[:, 0],
            name="Model",
            line=dict(color="black", dash="dot"),
            mode="lines",
        ),
    )
    fig = style_plot(
        fig, line_thickness=1.8, font_size=11, left_margin=20, buttom_margin=20
    )
    fig.write_image("plots/experiments/ct_mono_old.svg")


def plot_thiamine_gradient():
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
            range=[0, max(x)], showgrid=True, zeroline=True, dtick=12, title="Time [h]"
        ),
        yaxis=dict(
            range=[0, 0.4],
            showgrid=True,
            zeroline=True,
            dtick=0.1,
            title="OD",
        ),
        width=150,
        height=150,
        showlegend=False,
        # title="Thiamine gradient <i>O. anthropi</i>",
    )
    fig = style_plot(
        fig,
        line_thickness=0.8,
        font_size=11,
        left_margin=30,
        right_margin=0,
        buttom_margin=25,
        top_margin=0,
    )
    fig.write_image("plots/experiments/thiamine_gradient.svg")


def max_growth_rate_high_conc():
    df = pd.read_csv("/home/eric/ChiBioFlow/data/ct_phenotyping/metadata.csv")
    data = pd.read_csv("/home/eric/ChiBioFlow/data/ct_phenotyping/measurements.csv")
    slice = data[
        (data["240623_growth_phenotyping_ct_A1_time"] > 4.9)
        & (data["240623_growth_phenotyping_ct_A1_time"] < 9.99)
    ]
    slope, intercept, r_value, p_value, std_err = linregress(
        slice["240623_growth_phenotyping_ct_A1_time"],
        np.log(slice["240623_growth_phenotyping_ct_A1_measurement"]),
    )
    print("Slope Ct", slope)
    fig = go.Figure()
    lgs = df["linegroup"]
    for i, lg in enumerate(lgs):
        x = data[lg + "_time"]
        y = data[lg + "_measurement"]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="Ct",
                marker=dict(color=colors["ct"]),
                showlegend=(i == 0),
            )
        )
    fig.add_trace(
        go.Scatter(
            x=slice["240623_growth_phenotyping_ct_A1_time"],
            y=[
                np.exp(slope * i + intercept)
                for i in slice["240623_growth_phenotyping_ct_A1_time"]
            ],
            name="fit",
            line=dict(dash="dot", color="black"),
            mode="lines",
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)], showgrid=True, zeroline=True, dtick=8, title="Time [h]"
        ),
        yaxis=dict(
            range=[0, 1.3],
            showgrid=True,
            zeroline=True,
            dtick=0.2,
            title="OD",
        ),
        width=width,
        height=height,
    )
    fig = style_plot(fig, line_thickness=1.7, font_size=8)
    fig.write_image("plots/experiments/ct_high_conc.svg")

    fig = go.Figure()

    df = pd.read_csv("/home/eric/ChiBioFlow/data/at_oa/oa_high_conc/metadata.csv")

    data = pd.read_csv("/home/eric/ChiBioFlow/data/at_oa/oa_high_conc/measurements.csv")
    slice = data[
        (data["240623_growth_phenotyping_oa_A1_time"] > 4.9)
        & (data["240623_growth_phenotyping_oa_A1_time"] < 20.5)
    ]
    slope, intercept, r_value, p_value, std_err = linregress(
        slice["240623_growth_phenotyping_oa_A1_time"],
        np.log(slice["240623_growth_phenotyping_oa_A1_measurement"]),
    )
    print(slope)
    lgs = df["linegroup"]
    for i, lg in enumerate(lgs):
        x = data[lg + "_time"]
        y = data[lg + "_measurement"]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name="Oa",
                marker=dict(color=colors["oa"]),
                showlegend=(i == 0),
            )
        )
    fig.add_trace(
        go.Scatter(
            x=slice["240623_growth_phenotyping_oa_A1_time"],
            y=[
                np.exp(slope * i + intercept)
                for i in slice["240623_growth_phenotyping_oa_A1_time"]
            ],
            name="fit",
            mode="lines",
            line=dict(dash="dot", color="black"),
        )
    )
    fig.update_layout(
        xaxis=dict(
            range=[0, max(x)], showgrid=True, zeroline=True, dtick=8, title="Time [h]"
        ),
        yaxis=dict(
            range=[0, 1.3],
            showgrid=True,
            zeroline=True,
            dtick=0.2,
            title="OD",
        ),
        width=width,
        height=height,
    )
    fig = style_plot(
        fig,
        line_thickness=1.7,
        font_size=8,
        buttom_margin=20,
        top_margin=20,
        left_margin=20,
    )
    fig.write_image("plots/experiments/oa_high_conc.svg")


def spent_media_curves():
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/metadata.csv"
    )
    df_ct = df[df["species"] == "Comamonas testosteroni"]
    df_ct = df_ct.sort_values(by="exp_ID")
    df_oa = df[df["species"] == "Ochrobactrum anthropi"]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/measurements.csv"
    )
    fig = go.Figure()
    ys = []
    for j, cs in enumerate(["Spent media Ct", "Spent media Oa"]):
        line_counter = {"Batch 1": 0, "Batch 2": 0, "Batch 3": 0}
        for i, (lg, comment) in enumerate(
            df_ct[df_ct["carbon_source"] == cs][["linegroup", "comments"]].values
        ):
            if line_counter[comment] == 0:
                x = data[lg + "_time"]
                y = data[lg + "_measurement"]
                ys.append(y.to_numpy())
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        name=cs,
                        legendgroup=cs,
                        line=dict(color=colors["ct"]),
                        showlegend=(i == 0),
                        opacity=(0.5 if j == 0 else 1),
                    )
                )
                line_counter[comment] += 1
    ys = np.average(ys, axis=0)
    xs = x
    i, j = 0, 0
    for q, x in enumerate(xs):
        if x <= 24:
            i = q
        if x >= 60:
            j = q
            break
    slope, intercept, r_value, p_value, std_err = linregress(xs[i:j], np.log(ys[i:j]))
    print("Slope Ct", slope)
    """fig.add_trace(
        go.Scatter(
            x=xs[i:j],
            y=np.exp(slope * xs[i:j] + intercept),
            name="fit",
            line=dict(dash="dot", color=colors["ct"]),
        )
    )"""

    fig.update_layout(
        xaxis=dict(title="Time [h]", range=[0, 72], dtick=12),
        yaxis=dict(title="OD", range=[0, 0.08], dtick=0.02),
        showlegend=False,
        width=width,
        height=height,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        font_size=11,
        buttom_margin=25,
        left_margin=35,
        right_margin=0,
    )
    fig.write_image("plots/experiments/spent_media_ct.svg")
    fig = go.Figure()
    for j, cs in enumerate(["Spent media Ct", "Spent media Oa"]):
        line_counter = {"Batch 1": 0, "Batch 2": 0, "Batch 3": 0}
        for i, (lg, comment) in enumerate(
            df_oa[df_oa["carbon_source"] == cs][["linegroup", "comments"]].values
        ):
            if line_counter[comment] == 0:
                x = data[lg + "_time"]
                y = data[lg + "_measurement"]
                fig.add_trace(
                    go.Scatter(
                        x=x[5:],
                        y=y[5:],
                        name=cs,
                        legendgroup=cs,
                        mode="lines",
                        line=dict(
                            color=colors["oa"],
                        ),
                        opacity=(0.5 if j == 0 else 1),
                        showlegend=(i == 0),
                    )
                )
                line_counter[comment] += 1
    fig.update_layout(
        xaxis=dict(title="Time [h]", range=[0, 72], dtick=12),
        yaxis=dict(title="OD", range=[0, 0.03], dtick=0.01),
        title="<i>O. anthropi</i>",
        width=width,
        height=height,
        showlegend=False,
    )
    fig = style_plot(
        fig,
        line_thickness=0.8,
        font_size=8,
        buttom_margin=20,
        left_margin=25,
        top_margin=20,
        right_margin=0,
    )
    fig.write_image("plots/experiments/spent_media_oa.svg")


def Km_cufs():
    concentrations = {
        30: "#1f77b4",  # blue
        10: "#ff7f0e",  # orange
        1: "#2ca02c",  # green
        0.1: "#d62728",  # red
        0.01: "#9467bd",  # purple
        0.001: "#8c564b",  # brown
        0.0001: "#e377c2",  # pink
        "M9 + HMB": "#17becf",  # teal
    }

    slopes = {key: [] for key in list(concentrations.keys())}

    sheets = ["ct1", "ct2", "ct3", "oa1", "oa2", "oa3"]
    dfs = [
        pd.read_excel(
            "/home/eric/ChiBioFlow/data/at_oa/ct_oa_affinity_test/cfus.xlsx",
            sheet_name=sheet,
            index_col=0,
        )
        for sheet in sheets
    ]

    fig = go.Figure()

    sheets = ["ct1", "ct2", "ct3", "oa1", "oa2", "oa3"]
    dfs = [
        pd.read_excel(
            "/home/eric/ChiBioFlow/data/at_oa/ct_oa_affinity_test/cfus.xlsx",
            sheet_name=sheet,
            index_col=0,
        )
        for sheet in sheets
    ]

    x_values = list(dfs[0].columns)

    for c in list(concentrations.keys()):
        for i, df in enumerate(dfs[3:4]):  # first 3 replicates
            y_values = np.log(np.array(df.loc[c].values))
            slope, inter = linregress(x=range(len(y_values)), y=y_values)[:2]
            slopes[c].append(slope)
            fig.add_trace(
                go.Scatter(
                    x=x_values,
                    y=np.exp(y_values),
                    name=str(c),
                    marker=dict(color=concentrations[c]),
                    showlegend=(i == 0),
                    legendgroup=str(c),
                )
            )
            # Plot linear fit
            fig.add_trace(
                go.Scatter(
                    x=x_values,
                    y=np.exp([slope * x + inter for x in range(len(y_values))]),
                    line=dict(dash="dot", color=concentrations[c]),
                    showlegend=(i == 0),
                    legendgroup=str(c),
                )
            )
    fig.show()

    x = 1 / np.array(list(concentrations.keys())[3:-1])
    y = 1 / np.array([np.mean(slopes[c]) for c in list(concentrations.keys())[3:-1]])
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(color="black"),
        )
    )
    slope, inter = linregress(x, y)[:2]
    fig.add_trace(
        go.Scatter(
            x=x,
            y=[slope * i + inter for i in x],
            line=dict(dash="dot", color="black"),
        )
    )
    fig.write_image("tmp.svg")
    print(slopes)
    fig.show()
