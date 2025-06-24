from chibio_parser import (
    fluorescence_paresr,
    cfu_parser,
    ct_oa_parser,
    calibration_csv,
    calibration,
)
import pandas as pd
import numpy as np
from equations import *
import plotly.express as px


def get_cfus():
    replicats = ["replicate_1", "replicate_2", "replicate_3"]
    dfs = []
    sheets = []
    df = cfu_parser("/home/eric/ChiBioFlow/data/at_oa/241107_ct_oa")[0]
    df.insert(len(df.columns), "experiment", "ct_oa_thiamine")

    for i, r in enumerate(df["reactor"].unique()):
        sheet = pd.DataFrame(
            columns=["time", "CFUs/mL", "name", "species", "description", "figure"]
        )
        mask = df[(df["reactor"] == r) & (df["species"] == "ct")]
        (
            sheet["time"],
            sheet["CFUs/mL"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            mask["sample_time"],
            mask["average"],
            replicats[i],
            "Ct",
            "Ct Oa community chemostat experiment with 10000 nM thiamine",
            "1C",
        )
        sheets.append(sheet)
    for i, r in enumerate(df["reactor"].unique()):
        sheet = pd.DataFrame(
            columns=["time", "CFUs/mL", "name", "species", "description", "figure"]
        )
        mask = df[(df["reactor"] == r) & (df["species"] == "oa")]
        (
            sheet["time"],
            sheet["CFUs/mL"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            mask["sample_time"],
            mask["average"],
            replicats[i],
            "Oa",
            "Ct Oa community chemostat experiment with 10000 nM thiamine",
            "1C",
        )
        sheets.append(sheet)
    dfs.append(df)

    df = cfu_parser("/home/eric/ChiBioFlow/data/at_oa/241113_ct_oa")[0]
    df.insert(len(df.columns), "experiment", "ct_oa")
    dfs.append(df)
    for i, r in enumerate(df["reactor"].unique()):
        sheet = pd.DataFrame(
            columns=["time", "CFUs/mL", "name", "species", "description", "figure"]
        )
        mask = df[(df["reactor"] == r) & (df["species"] == "ct")]
        (
            sheet["time"],
            sheet["CFUs/mL"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            mask["sample_time"],
            mask["average"],
            replicats[i],
            "Ct",
            "Ct Oa cross-feeding community chemostat experiment",
            "1D",
        )
        sheets.append(sheet)
    for i, r in enumerate(df["reactor"].unique()):
        sheet = pd.DataFrame(
            columns=["time", "CFUs/mL", "name", "species", "description", "figure"]
        )
        mask = df[(df["reactor"] == r) & (df["species"] == "oa")]
        (
            sheet["time"],
            sheet["CFUs/mL"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            mask["sample_time"],
            mask["average"],
            replicats[i],
            "Oa",
            "Ct Oa cross-feeding community chemostat experiment",
            "1D",
        )
        sheets.append(sheet)
    pd.concat(sheets).to_excel(
        "../data/data.xlsx", index=False, sheet_name="Chemostat CFU data"
    )
    return pd.concat(dfs)


def get_od_chemostats():
    dfs = []
    Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/250310_oa_mono")
    df = calibration_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250310_oa_mono/calibration.csv", Ms
    )
    df.insert(len(df.columns), "experiment", "oa_mono")
    dfs.append(df)
    Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/250327_oa_mono")
    df = calibration_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250327_oa_mono/calibration.csv", Ms
    )
    df.insert(len(df.columns), "experiment", "oa_mono_repeat")
    dfs.append(df)
    Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/250320_ct_mono")
    df = calibration_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250320_ct_mono/calibration.csv", Ms
    )
    df.insert(len(df.columns), "experiment", "ct_mono")
    dfs.append(df)
    Ms = fluorescence_paresr("/home/eric/ChiBioFlow/data/at_oa/241101_ct_mono")
    Ms = [Ms[Ms["reactor"] == r] for r in ["M0", "M1", "M2"]]
    Ms = [M[M["exp_time"] < 500] for M in Ms]
    k_eq, b_eq = calibration()
    for M in Ms[1:]:
        M.index = range(len(M))
        R0_raw = np.average(M.loc[3:5]["od_measured"])
        R1_raw = np.average(M.loc[22.5 * 60 : 23.5 * 60]["od_measured"])
        k_num = float(k_eq.subs({"OD1": 0.05, "OD2": 0.141, R1: R0_raw, R2: R1_raw}))
        b_num = float(b_eq.subs({"OD1": 0.05, "OD2": 0.141, R1: R0_raw, R2: R1_raw}))
        OD = k_num * np.log10(M["od_measured"]) + b_num
        M["od_calibrated"] = OD
        M["experiment"] = "ct_mono_old"
        dfs.append(M)
    return pd.concat(dfs)


def ct_oa_plate_reader():
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250328_ct_oa_thiamine_gradient/data/metadata.csv"
    )
    df_ct = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
        & (df["species"] == "Comamonas testosteroni")
        & (df["comments"] == "10000 nM thiamine")
    ]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250328_ct_oa_thiamine_gradient/data/measurements.csv"
    )
    sheets = []

    replicates = ["replicate_1", "replicate_2", "replicate_3"]
    for i, lg in enumerate(df_ct["linegroup"]):
        sheet = pd.DataFrame(
            columns=["time", "OD", "name", "species", "description", "figure"]
        )
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        (
            sheet["time"],
            sheet["OD"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            x,
            y,
            replicates[i],
            "Ct",
            "Ct mono-culture growth curve with 10000 nM thiamine",
            "1A",
        )
        sheets.append(sheet)

    df_oa = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
        & (df["comments"] == "10000 nM thiamine")
    ]

    for i, lg in enumerate(df_oa["linegroup"]):
        sheet = pd.DataFrame(
            columns=["time", "OD", "name", "species", "description", "figure"]
        )
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        (
            sheet["time"],
            sheet["OD"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            x,
            y,
            replicates[i],
            "Oa",
            "Oa mono-culture growth curve with 10000 nM thiamine",
            "1A",
        )
        sheets.append(sheet)

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
        sheet = pd.DataFrame(
            columns=["time", "OD", "name", "species", "description", "figure"]
        )
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        (
            sheet["time"],
            sheet["OD"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            x,
            y,
            replicates[i],
            "Oa",
            "Oa mono-culture growth curve with no thiamine",
            "1A",
        )
        sheets.append(sheet)

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
        sheet = pd.DataFrame(
            columns=["time", "OD", "name", "species", "description", "figure"]
        )
        x = data[lg + "_time"][data[lg + "_time"] < 36]
        y = data[lg + "_measurement"][: len(x)]
        (
            sheet["time"],
            sheet["OD"],
            sheet["name"],
            sheet["species"],
            sheet["description"],
            sheet["figure"],
        ) = (
            x,
            y,
            replicates[i],
            "Oa",
            "Oa mono-culture growth curve in spent-media of Ct",
            "1A",
        )
        sheets.append(sheet)

        out = pd.concat(sheets)
    fig = px.line(
        out,
        x="time",
        y="OD",
        color="species",
        facet_col="description",
        line_group="name",
    )

    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/metadata.csv"
    )
    df = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_oa_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
    ]
    mask = []
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
    for c in df["comments"]:
        if c in keep:
            mask.append(True)
        else:
            mask.append(False)
    df = df[mask]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250313_oa_thiamine_gradient/data/measurements.csv"
    )

    for t_conc in df["comments"].unique():
        mask = df["comments"] == t_conc
        for i, lg in enumerate(df[mask]["linegroup"]):
            sheet = pd.DataFrame(
                columns=["time", "OD", "name", "species", "description", "figure"]
            )
            x = data[lg + "_time"][data[lg + "_time"] < 36]
            y = data[lg + "_measurement"][: len(x)]
            (
                sheet["time"],
                sheet["OD"],
                sheet["name"],
                sheet["species"],
                sheet["description"],
                sheet["figure"],
            ) = (
                x,
                y,
                replicates[i],
                "Oa",
                "Oa mono-culture growth curve with " + t_conc,
                "1A",
            )
            sheets.append(sheet)
    df = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/metadata.csv"
    )
    df_ct = df[df["species"] == "Comamonas testosteroni"]
    df_ct = df_ct.sort_values(by="exp_ID")
    df_oa = df[df["species"] == "Ochrobactrum anthropi"]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250411_spent_media_growth_curves/measurements.csv"
    )
    names = [
        "Ct grown in Ct's spent chemostat media",
        "Ct grown in Oa's spent chemostat media",
    ]

    for j, cs in enumerate(["Spent media Ct", "Spent media Oa"]):

        line_counter = {"Batch 1": 0, "Batch 2": 0, "Batch 3": 0}
        for i, (lg, comment) in enumerate(
            df_ct[df_ct["carbon_source"] == cs][["linegroup", "comments"]].values
        ):
            sheet = pd.DataFrame(
                columns=["time", "OD", "name", "species", "description", "figure"]
            )
            if line_counter[comment] == 0:
                x = data[lg + "_time"]
                y = data[lg + "_measurement"]
                (
                    sheet["time"],
                    sheet["OD"],
                    sheet["name"],
                    sheet["species"],
                    sheet["description"],
                    sheet["figure"],
                ) = (
                    x,
                    y,
                    "replicate_" + comment.split(" ")[-1],
                    "Ct",
                    names[j],
                    "2C",
                )
                line_counter[comment] += 1

                sheets.append(sheet)
    names = [
        "Oa grown in Oa's spent chemostat media",
        "Oa grown in Ct's spent chemostat media",
    ]
    for j, cs in enumerate(["Spent media Ct", "Spent media Oa"]):
        line_counter = {"Batch 1": 0, "Batch 2": 0, "Batch 3": 0}
        for i, (lg, comment) in enumerate(
            df_oa[df_oa["carbon_source"] == cs][["linegroup", "comments"]].values
        ):
            sheet = pd.DataFrame(
                columns=["time", "OD", "name", "species", "description", "figure"]
            )
            if line_counter[comment] == 0:

                x = data[lg + "_time"]
                y = data[lg + "_measurement"]
                (
                    sheet["time"],
                    sheet["OD"],
                    sheet["name"],
                    sheet["species"],
                    sheet["description"],
                    sheet["figure"],
                ) = (
                    x,
                    y,
                    "replicate_" + comment.split(" ")[-1],
                    "Oa",
                    names[j],
                    None,
                )
                line_counter[comment] += 1

                sheets.append(sheet)

    pd.concat(sheets).to_excel(
        "../data/data.xlsx", index=False, sheet_name="Monoculture batch growth curves"
    )
    return pd.concat(sheets)
