from chibio_parser import fluorescence_paresr
from chibio_parser import calibration_csv
import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *


def read_raw(plot=False):
    raw = fluorescence_paresr("260629_oa_washout_ct_oa_no_cs")
    raw = raw[(raw["exp_time"] >= 0.1) & (raw["exp_time"] <= 93.5)]
    if plot:
        fig = px.line(raw, x="exp_time", y="od_measured", color="reactor")
        fig.show()
    return raw


def write_cal():
    raw = read_raw()
    od0_plate_reader = {"M0": 0.27, "M1": 0.14, "M2": 0.15}
    od1_plate_reader = {"M0": 0.013, "M1": 0.02, "M2": 0}

    cal = pd.DataFrame(
        columns=["reactor", "OD0", "R0", "OD1", "R1"],
        index=["M0", "M1", "M2"],
    )
    cal["reactor"] = cal.index

    for reactor_name, df in raw.groupby("reactor"):
        R0 = np.average(df[df["exp_time"] <= 0.2]["od_measured"])
        R1 = np.average(df[df["exp_time"] >= 92]["od_measured"])
        row = [
            reactor_name,
            od0_plate_reader[reactor_name],
            R0,
            od1_plate_reader[reactor_name],
            R1,
        ]
        cal.loc[reactor_name] = row
    cal.to_csv("calibration.csv", index=False)


# write_cal()
df = calibration_csv(
    "/home/eric/ChiBioFlow/data/260629_oa_washout_ct_oa_no_cs/calibration.csv",
    read_raw(),
)
species = ["oa", "oa", "ct"]
titles = ["Oa, no thiamine", "Oa, no cs", "Ct, no cs"]
fnames = ["oa_no_thiamine.svg", "oa_no_cs.svg", "ct_no_cs.svg"]
ods_engel = pd.read_excel("ods/engel/ods.xlsx")
for i, (reactor_name, r_df) in enumerate(df.groupby("reactor")):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=r_df["exp_time"][::15],
            y=r_df["od_calibrated"][::15],
            marker=dict(color=colors[species[i]]),
            name="Chi.Bio",
            opacity=0.6,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=ods_engel["time"],
            y=ods_engel[reactor_name],
            mode="markers",
            marker=dict(color=colors[species[i]]),
            name="platereader",
        )
    )

    fig.update_layout(
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD", ticks="inside"),
        title=titles[i],
        width=250,
        height=200,
        legend=dict(
            xref="paper", yref="paper", xanchor="right", yanchor="top", x=0.99, y=0.99
        ),
    )
    fig = style_plot(fig, marker_size=6, line_thickness=2, top_margin=30)
    fig.write_image(fnames[i])
