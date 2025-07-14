import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from style import *
import math
import scipy.stats as stats


media = ["SM_042025_MPTA_b1_1", "SM_042025_MPTA_b2_2", "SM_042025_MPTA_b3_3"]
ct_c = [
    "SM_042025_MPTA_ct_c_b2_5",
    "SM_042025_MPTA_ct_c_b3_6",
    "SM_042025_MPTA_ct_c_b1_4",
]
oa_c = [
    "SM_042025_MPTA_oa_c_b3_9",
    "SM_042025_MPTA_oa_c_b1_7",
    "SM_042025_MPTA_oa_c_b2_8",
]
ct_b = [
    "SM_042025_MPTA_ct_b_b2_11",
    "SM_042025_MPTA_ct_b_b3_12",
    "SM_042025_MPTA_ct_b_b1_10",
]
oa_b = [
    "SM_042025_MPTA_oa_b_b1_13",
    "SM_042025_MPTA_oa_b_b2_14",
    "SM_042025_MPTA_oa_b_b3_15",
]

raw = pd.read_excel(
    "../data/250610_ms_data/raw_data_absolut.xlsx",
)
raw.index = raw["metabolite"]
raw = raw.replace("< LOQ", None)

meta = pd.read_excel("../data/250610_ms_data/meta.xlsx")
meta.index = meta["metabolite"]

raw.insert(len(raw.columns), "group", None)
for m in raw["metabolite"]:
    raw.loc[m, "group"] = meta.loc[m]["group"]

symbols = {"Ct": "cross", "Oa": "square"}


def plotting():
    for group in set(raw["group"]):
        df = raw[raw["group"] == group]
        f = group.replace(" ", "_").replace("/", "_") + ".pdf"
        rows = math.ceil(len(df) / 4)
        fig = make_subplots(rows=rows, cols=4, subplot_titles=df["metabolite"].values)

        for i, m in enumerate(df["metabolite"]):
            row = i // 4 + 1
            col = i % 4 + 1

            x = 3 * ["Media"]
            y = df.loc[m][media].values

            fig.add_trace(
                go.Box(
                    x=x,
                    y=y,
                    marker=dict(color="black"),
                    boxpoints="all",
                    jitter=1,
                    pointpos=0,
                    boxmean=True,
                    quartilemethod="linear",
                ),
                row=row,
                col=col,
            )

            x = 3 * ["Ct Chemostat"]
            y = df.loc[m][ct_c].values
            fig.add_trace(
                go.Box(
                    x=x,
                    y=y,
                    marker=dict(color=colors["ct"]),
                    boxpoints="all",
                    jitter=1,
                    pointpos=0,
                    boxmean=True,
                    quartilemethod="linear",
                ),
                row=row,
                col=col,
            )

            x = 3 * ["Oa Chemostat"]
            y = df.loc[m][oa_c].values
            fig.add_trace(
                go.Box(
                    x=x,
                    y=y,
                    marker=dict(color=colors["oa"]),
                    boxpoints="all",
                    jitter=1,
                    pointpos=0,
                    boxmean=True,
                    quartilemethod="linear",
                ),
                row=row,
                col=col,
            )
            x = 3 * ["Ct in Oa"]
            y = df.loc[m][ct_b].values
            fig.add_trace(
                go.Box(
                    x=x,
                    y=y,
                    marker=dict(color=colors["ct"]),
                    boxpoints="all",
                    jitter=1,
                    pointpos=0,
                    boxmean=True,
                    quartilemethod="linear",
                ),
                row=row,
                col=col,
            )
            x = 3 * ["Oa in Ct"]
            y = df.loc[m][oa_b].values
            fig.add_trace(
                go.Box(
                    x=x,
                    y=y,
                    marker=dict(color=colors["oa"]),
                    boxpoints="all",
                    jitter=1,
                    pointpos=0,
                    boxmean=True,
                    quartilemethod="linear",
                ),
                row=row,
                col=col,
            )
        fig.update_layout(showlegend=False, width=800, height=rows * 150),
        fig.update_yaxes(type="log", nticks=4, title=df.loc[m]["unit"])
        fig = style_plot(fig, line_thickness=1, marker_size=5, font_size=11)
        fig.write_image("plots/ms_analysis/absolut/" + f)


def media_metabolites():
    relevant = ["Isoleucine", "Glutamine", "Lactate", "Citrate"]
    fig = go.Figure()
    raw = raw[["metabolite", "group"] + media].copy()
    for m, row in raw.iterrows():
        if m not in relevant:
            continue
        fig.add_trace(
            go.Box(
                x=[m, m, m],
                y=row[media].values,
                name=m,
                marker=dict(color="black"),
                boxpoints="all",
                jitter=1,
                pointpos=0,
                boxmean=True,
                quartilemethod="linear",
                showlegend=False,
            ),
        )
    fig.update_yaxes(type="log")
    fig.update_layout(
        xaxis_title="Metabolite",
        yaxis_title="Concentration [µM]",
        width=width,
        height=1.5 * height,
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        marker_size=5,
        font_size=11,
        buttom_margin=80,
        top_margin=0,
        left_margin=40,
        right_margin=0,
    )
    fig.write_image("plots/ms_analysis/absolut/media.svg")


def fig2f():

    df = raw[raw["metabolite"] == "Cis-Aconitate"]
    fig = go.Figure()
    fig.add_trace(
        go.Box(
            x=["Media"],
            y=[],
            name="Media",
            text="LOQ",
            marker=dict(color="black"),
            boxpoints="all",
            jitter=1,
            pointpos=0,
            boxmean=True,
            quartilemethod="linear",
            showlegend=False,
        ),
    )
    fig.add_annotation(
        x="Media",
        y=0,
        text="< LOQ",
        showarrow=False,
        textangle=90,  #
        # font=dict(size=14, color="gray"),
        align="center",
        bgcolor="white",
        # bordercolor="black",
        borderwidth=1,
        yshift=-30,  # shift the text box above baseline
    )
    fig.add_trace(
        go.Box(
            x=3 * ["Ct Chemostat"],
            y=df[ct_c].values[0],
            name="Ct Chemostat",
            marker=dict(color=colors["ct"]),
            boxpoints="all",
            jitter=1,
            pointpos=0,
            boxmean=True,
            quartilemethod="linear",
            showlegend=False,
        ),
    )
    fig.add_trace(
        go.Box(
            x=3 * ["Oa Chemostat"],
            y=df[oa_c].values[0],
            name="Oa Chemostat",
            marker=dict(color=colors["oa"]),
            boxpoints="all",
            jitter=1,
            pointpos=0,
            boxmean=True,
            quartilemethod="linear",
            showlegend=False,
        ),
    )
    fig.add_trace(
        go.Box(
            x=3 * ["Ct in Oa"],
            y=[],
            name="Ct in Oa",
            marker=dict(color=colors["ct"]),
            boxpoints="all",
            jitter=1,
            pointpos=0,
            boxmean=True,
            quartilemethod="linear",
            showlegend=False,
        ),
    )
    fig.add_annotation(
        x="Ct in Oa",
        y=0,
        text="< LOQ",
        showarrow=False,
        textangle=90,  #
        # font=dict(size=14, color="gray"),
        align="center",
        bgcolor="white",
        # bordercolor="black",
        borderwidth=1,
        yshift=-30,  # shift the text box above baseline
    )

    fig.add_trace(
        go.Box(
            x=3 * ["Oa in Ct"],
            y=df[oa_b].values[0],
            name="Oa in Ct",
            marker=dict(color=colors["oa"]),
            boxpoints="all",
            jitter=1,
            pointpos=0,
            boxmean=True,
            quartilemethod="linear",
            showlegend=False,
        ),
    )

    fig.update_yaxes(title="Concentration [µM]", type="log")
    fig.update_layout(
        yaxis_title="Concentration [µM]",
        width=width * 0.8,
        height=1.3 * height,
        title="Cis-Aconitate",
    )
    fig = style_plot(
        fig,
        line_thickness=1,
        marker_size=5,
        font_size=11,
        buttom_margin=30,
        top_margin=30,
        left_margin=30,
        right_margin=0,
    )
    fig.write_image("plots/ms_analysis/absolut/fig2f.svg")
