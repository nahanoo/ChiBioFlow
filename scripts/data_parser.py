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
    dfs = []
    df = cfu_parser("/home/eric/ChiBioFlow/data/at_oa/241107_ct_oa")[0]
    df.insert(len(df.columns), "experiment", "ct_oa_thiamine")
    dfs.append(df)

    df = cfu_parser("/home/eric/ChiBioFlow/data/at_oa/241113_ct_oa")[0]
    df.insert(len(df.columns), "experiment", "ct_oa")
    dfs.append(df)
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
    df_oa = df[
        (df["exp_ID"] == "ct_oa_chemostat_project/_thiamine_gradient")
        & (df["species"] == "Ochrobactrum anthropi")
        & (df["comments"] == "10000 nM thiamine")
    ]
    data = pd.read_csv(
        "/home/eric/ChiBioFlow/data/at_oa/250328_ct_oa_thiamine_gradient/data/measurements.csv"
    )

    dfs = []

    metas = [df_ct, df_oa]
    for meta in metas:
        for i, lg in enumerate(meta["linegroup"]):
            time = data[lg + "_time"]
            od = data[lg + "_measurement"]
            df = pd.DataFrame(columns=["time", "species", "linegroup", "od"])
            df["time"], df["species"], df["linegroup"], df["od"] = time, "ct", lg, od
            dfs.append(df)
    out = pd.concat(dfs)
    return out[out["species"] == "ct"], out[out["species"] == "oa"]
