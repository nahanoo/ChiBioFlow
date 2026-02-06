import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from style import *
import math
import scipy.stats as stats
from plotly.subplots import make_subplots


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

raw = pd.read_excel("../data/250610_ms_data/raw_data.xlsx")
raw.index = raw["metabolite"]

meta = pd.read_excel("../data/250610_ms_data/meta.xlsx")
meta.index = meta["metabolite"]

raw.insert(len(raw.columns), "group", None)
for m in raw["metabolite"]:
    raw.loc[m, "group"] = meta.loc[m]["group"]

colors = {
    group: color
    for group, color in zip(
        set(meta["group"]),
        [
            "#1f77b4",
            "#ff7f0e",
            "#e377c2",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#2ca02c",
            "#7f7f7f",
        ],
    )
}
symbols = {"Ct": "cross", "Oa": "square"}


def plot_metabolic_classes():
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
        fig.update_yaxes(type="log", nticks=4)
        fig = style_plot(fig, line_thickness=1, marker_size=5, font_size=11)
        fig.write_image("plots/ms_analysis/relativ/" + f)


pairs = {
    "media_ct_c": [media, ct_c],
    "media_oa_c": [media, oa_c],
    "ct_c_oa_b": [ct_c, oa_b],
    "oa_c_ct_b": [oa_c, ct_b],
}
df = pd.DataFrame(
    columns=[
        "metabolite",
        "group",
        "media_ct_c",
        "media_oa_c",
        "ct_c_oa_b",
        "oa_c_ct_b",
    ]
)
for i, m in enumerate(raw["metabolite"]):
    for key, pair in pairs.items():
        group = meta.loc[m]["group"]
        fold_change = np.log2(
            np.mean(raw.loc[m][pair[1]].values) / np.mean(raw.loc[m][pair[0]].values)
        )
        df.at[i, key] = [
            stats.ttest_ind(
                list(raw.loc[m][pair[1]].values),
                list(raw.loc[m][pair[0]].values),
                equal_var=False,
            )[:2],
            fold_change,
        ]

        df.at[i, "metabolite"] = m
        df.at[i, "group"] = group
df.index = df["metabolite"]
df = df.sort_values(by="group")


def leakage_consumption():
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=["Ct", "Oa", "Ct", "Oa"],
        vertical_spacing=0.12,
        shared_yaxes=True,
        shared_xaxes=True,
    )
    groups = []
    consumed = []
    for i, (m, ct_batch) in enumerate(zip(df["metabolite"], df["oa_c_ct_b"])):
        if (ct_batch[0][1] < 0.05) & (ct_batch[1] < 0):
            fig.add_trace(
                go.Scatter(
                    x=[ct_batch[1]],
                    y=-np.log10([ct_batch[0][1]]),
                    marker=dict(
                        color=colors_metabolites[meta.loc[m]["group"]],
                    ),
                    hovertext=[meta.loc[m]["group"] + "<br>" + m + "<br>Ct"],
                    textfont=dict(size=8),
                    textposition="middle right",
                    mode="markers+text",
                    text=m,
                    showlegend=False,
                ),
                row=2,
                col=1,
            )
            groups.append(meta.loc[m]["group"])
            consumed.append(m)

    for i, (m, oa_batch) in enumerate(zip(df["metabolite"], df["ct_c_oa_b"])):
        if (oa_batch[0][1] < 0.05) & (oa_batch[1] < 0):
            fig.add_trace(
                go.Scatter(
                    x=[oa_batch[1]],
                    y=-np.log10([oa_batch[0][1]]),
                    marker=dict(
                        color=colors_metabolites[meta.loc[m]["group"]],
                    ),
                    hovertext=[meta.loc[m]["group"] + "<br>" + m],
                    showlegend=False,
                    textfont=dict(size=8),
                    textposition="middle right",
                    mode="markers+text",
                    text=m,
                ),
                row=2,
                col=2,
            )
            groups.append(meta.loc[m]["group"])
            consumed.append(m)

    for group in sorted(list(set(groups))):
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(color=colors_metabolites[group]),
                name=group,
            )
        )
    for i, (m, media_ct) in enumerate(zip(df["metabolite"], df["media_ct_c"])):
        if media_ct[0][1] < 0.05:
            fig.add_trace(
                go.Scatter(
                    x=[media_ct[1]],
                    y=-np.log10([media_ct[0][1]]),
                    marker=dict(
                        color=colors_metabolites[meta.loc[m]["group"]],
                    ),
                    hovertext=[meta.loc[m]["group"] + "<br>" + m + "<br>Ct"],
                    showlegend=False,
                    textfont=dict(size=8),
                    textposition="middle left",
                    mode=("markers+text" if m in consumed else "markers"),
                    text=m,
                ),
                row=1,
                col=1,
            )

    for i, (m, media_oa) in enumerate(zip(df["metabolite"], df["media_oa_c"])):
        if media_oa[0][1] < 0.05:
            fig.add_trace(
                go.Scatter(
                    x=[media_oa[1]],
                    y=-np.log10([media_oa[0][1]]),
                    marker=dict(
                        color=colors_metabolites[meta.loc[m]["group"]],
                    ),
                    hovertext=[meta.loc[m]["group"] + "<br>" + m],
                    showlegend=False,
                    textfont=dict(size=8),
                    textposition="middle left",
                    mode=("markers+text" if m in consumed else "markers"),
                    text=m,
                ),
                row=1,
                col=2,
            )
            groups.append(meta.loc[m]["group"])

    fig.update_layout(
        width=3 * width,
        height=2.6 * height,
    )
    fig.for_each_xaxis(lambda x: x.update(ticks="inside"))
    fig.for_each_yaxis(lambda y: y.update(ticks="inside"))
    fig = style_plot(
        fig,
        line_thickness=1,
        marker_size=5,
        font_size=11,
        buttom_margin=30,
        top_margin=40,
        left_margin=30,
        right_margin=30,
    )
    fig["layout"]["xaxis3"]["title"] = "log<sub>2</sub> fold change"
    fig["layout"]["yaxis1"]["title"] = "-log<sub>10</sub> p-value"
    # fig.show()
    fig.write_image("plots/ms_analysis/relativ/fig3a.svg")


leakage_consumption()


def fig2e():
    pass


def sfig2b(raw):
    fig = go.Figure()
    raw = raw[["metabolite", "group"] + media]
    df = raw[media][raw[media] > 10000].dropna()
    for m, row in df.iterrows():
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
        yaxis_title="Peak area",
        width=2.5 * width,
        height=1.5 * height,
    )
    fig = style_plot(
        fig,
        marker_size=5,
        font_size=11,
        buttom_margin=10,
        top_margin=0,
        left_margin=40,
        right_margin=30,
    )
    fig.write_image("plots/ms_analysis/relativ/sfig2b.svg")


def sfig2c(raw):
    colors = {
        "blue": "#000080",
        "ct": "#7570B3",
        "oa": "#D95F02",
        "ms": "#E6AB02",
        "at": "#1B9E77",
        "Spent media Ct": "#1B9E77",
        "Spent media Oa": "#E7298A",
        "H20": "gray",
    }

    t = raw[raw["metabolite"] == "Thiamine"]
    fig = go.Figure()
    for i, m in enumerate(t["metabolite"]):
        x = 3 * ["Media"]
        y = t.loc[m][media].values

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
        )

        x = 3 * ["Ct Chemostat"]
        y = t.loc[m][ct_c].values
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
        )

        x = 3 * ["Oa Chemostat"]
        y = t.loc[m][oa_c].values
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
        )
        x = 3 * ["Ct in Oa"]
        y = t.loc[m][ct_b].values
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
        )
        x = 3 * ["Oa in Ct"]
        y = t.loc[m][oa_b].values
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
        )
    fig.update_layout(
        showlegend=False, width=width, height=height * 1.5, yaxis=dict(type="log")
    ),
    fig.update_yaxes(type="log", nticks=4)
    fig = style_plot(fig, marker_size=5, font_size=11)
    fig.write_image("plots/ms_analysis/relativ/sfig2c.svg")

    fig = go.Figure()
    for i, m in enumerate(t["metabolite"]):
        x = 3 * ["Media"]
        y = t.loc[m][media].values
        y_media = y
        fig.add_trace(
            go.Box(
                x=x,
                y=10 / y_media * y,
                marker=dict(color="black"),
                boxpoints="all",
                jitter=1,
                pointpos=0,
                boxmean=True,
                quartilemethod="linear",
            ),
        )

        x = 3 * ["Ct Chemostat"]
        y = t.loc[m][ct_c].values
        fig.add_trace(
            go.Box(
                x=x,
                y=10 / y_media * y,
                marker=dict(color=colors["ct"]),
                boxpoints="all",
                jitter=1,
                pointpos=0,
                boxmean=True,
                quartilemethod="linear",
            ),
        )

        x = 3 * ["Oa Chemostat"]
        y = t.loc[m][oa_c].values
        fig.add_trace(
            go.Box(
                x=x,
                y=10 / y_media * y,
                marker=dict(color=colors["oa"]),
                boxpoints="all",
                jitter=1,
                pointpos=0,
                boxmean=True,
                quartilemethod="linear",
            ),
        )
        x = 3 * ["Ct in Oa"]
        y = t.loc[m][ct_b].values
        fig.add_trace(
            go.Box(
                x=x,
                y=10 / y_media * y,
                marker=dict(color=colors["ct"]),
                boxpoints="all",
                jitter=1,
                pointpos=0,
                boxmean=True,
                quartilemethod="linear",
            ),
        )
        x = 3 * ["Oa in Ct"]
        y = t.loc[m][oa_b].values
        fig.add_trace(
            go.Box(
                x=x,
                y=10 / y_media * y,
                marker=dict(color=colors["oa"]),
                boxpoints="all",
                jitter=1,
                pointpos=0,
                boxmean=True,
                quartilemethod="linear",
            ),
        )
    fig.update_layout(showlegend=False, width=width, height=height),
    fig.update_yaxes(type="log", nticks=4)
    fig = style_plot(fig, line_thickness=1, marker_size=5, font_size=11)
    fig.write_image("plots/ms_analysis/relativ/thiamine_absolut.svg")


sfig2c(raw)
