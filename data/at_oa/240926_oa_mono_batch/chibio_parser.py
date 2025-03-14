from statistics import stdev, mean
from os.path import join
import pandas as pd
from glob import glob
import json
import numpy as np


def cfu_parser(e):
    with open(
        join("/", "home", "eric", "ChiBioFlow", "data", e, "order.txt"), "r"
    ) as f:
        order = f.read().rstrip().split(",")

    """Parses counted CFUs based on template xlsx"""
    # For eveery species there is one xlsx sheet
    fs = glob(join("/", "home", "eric", "ChiBioFlow", "data", e, "cfus*.xlsx"))
    # df for concatenating species sheets
    dfs = []
    for f in fs:
        # Species name
        species = f[-7:-5]
        cfus = pd.read_excel(f, header=1)
        cfus.insert(0, "species", species)
        dfs.append(cfus)
    # Reindexing because of concat
    df = pd.concat(dfs)
    df.index = range(len(df))

    # Splitting "|" separated triplicates
    counts = []
    for count, dilution in zip(df["count"], df["dilution"]):
        # Conversion to CFUs/mL sample volume 5 uL
        counts.append([(int(n) / 10**dilution) * 200 for n in count.split("|")])
    # Calculating average and stdev
    df["average"] = [mean(count) for count in counts]
    df["stdev"] = [stdev(count) for count in counts]

    # Adding total counts for composition
    df.insert(len(df.columns), "total", None)
    for t in df["sample_time"]:
        # Subsetting for reactor and sample time
        for reactor in order:
            tmp = df[df["sample_time"] == t]
            tmp = tmp[tmp["reactor"] == reactor]
            total = tmp["average"].sum()
            for i in tmp.index:
                df.at[i, "total"] = total
    # Adding composition
    df.insert(len(df.columns), "composition", None)
    df["composition"] = df["average"] / df["total"]
    return df, order


def chibio_parser(e, down_sample=False, sample_size=10):
    def average(df):
        out = pd.DataFrame(
            columns=["exp_time", "od_measured", "FP1_base", "FP1_emit1", "FP1_emit2"]
        )
        i = 0
        j = 0
        while True:
            try:
                sliece = df.loc[i : i + sample_size]
                out.loc[j] = [
                    sliece.iloc[-1]["exp_time"],
                    np.average(sliece["od_measured"]),
                    np.average(sliece["FP1_base"]),
                    np.average(sliece["FP1_emit1"]),
                    np.average(sliece["FP1_emit2"]),
                ]
                i += sample_size
                j += 1
            except IndexError:
                break
        return out

    with open(
        join("/", "home", "eric", "ChiBioFlow", "data", e, "order.txt"), "r"
    ) as f:
        order = f.read().rstrip().split(",")
    dfs = []
    for reactor in order:
        f = glob(
            join("/", "home", "eric", "ChiBioFlow", "data", e, reactor, "*data.csv")
        )[0]
        data = pd.read_csv(
            f, usecols=["exp_time", "od_measured", "FP1_base", "FP1_emit1", "FP1_emit2"]
        )
        if down_sample:
            data = average(data)
        data.insert(1, "reactor", reactor)
        data["exp_time"] = data["exp_time"] / 60 / 60
        dfs.append(data)

    df = pd.concat(dfs)
    dfs = []
    for c in ["od_measured", "FP1_emit1", "FP1_emit2", "FP1_base"]:
        tmp = pd.DataFrame(columns=["exp_time", "reactor", "measurement", "sensor"])
        tmp["reactor"] = df["reactor"]
        tmp["exp_time"] = df["exp_time"]
        tmp["measurement"] = df[c]
        tmp["sensor"] = c
        dfs.append(tmp)
    df = pd.concat(dfs)

    return df, order


def fluorescence_paresr(
    e,
    od_led_channel=False,
    od_led_converter=False,
    down_sample=True,
    window_size=10,
    od_filter=None,
):
    with open(
        join("/", "home", "eric", "ChiBioFlow", "data", e, "order.txt"), "r"
    ) as f:
        order = f.read().rstrip().split(",")
    dfs = []
    for reactor in order:
        f = glob(
            join("/", "home", "eric", "ChiBioFlow", "data", e, reactor, "*data.csv")
        )[0]
        columns = [
            "exp_time",
            "od_measured",
            "FP1_base",
            "FP1_emit1",
            "FP1_emit2",
            "FP2_base",
            "FP2_emit1",
            "FP2_emit2",
            "FP3_base",
            "FP3_emit1",
            "FP3_emit2",
            "reactor",
            "raw_FP1_emit1",
            "raw_FP1_emit2",
            "raw_FP2_emit1",
            "raw_FP2_emit2",
            "raw_FP3_emit1",
            "raw_FP3_emit2",
            "od_measured_led",
        ]
        data_columns = columns[:11]
        data = pd.read_csv(
            f,
            usecols=[
                "exp_time",
                "od_measured",
                "FP1_base",
                "FP1_emit1",
                "FP1_emit2",
                "FP2_base",
                "FP2_emit1",
                "FP2_emit2",
                "FP3_base",
                "FP3_emit1",
                "FP3_emit2",
            ],
        )
        df = pd.DataFrame(columns=columns)
        for c in data_columns:
            df[c] = data[c]
        FP1_emit = columns[3:5]
        FP1_emit_raw = columns[12:14]
        for emit, raw in zip(FP1_emit, FP1_emit_raw):
            df[raw] = data[emit] * data["FP1_base"]

        FP2_emit = columns[6:8]
        FP2_emit_raw = columns[14:16]
        for emit, raw in zip(FP2_emit, FP2_emit_raw):
            df[raw] = data[emit] * data["FP2_base"]

        FP3_emit = columns[9:11]
        FP3_emit_raw = columns[16:18]
        for emit, raw in zip(FP3_emit, FP3_emit_raw):
            df[raw] = data[emit] * data["FP3_base"]

        df["reactor"] = reactor
        df["exp_time"] = df["exp_time"] / 60 / 60

        if od_led_converter:
            df["od_measured_led"] = data[od_led_channel] * od_led_converter

        dfs.append(df)

    out = pd.concat(dfs)

    if down_sample:
        dfs = []
        for reactor in set(out["reactor"]):
            r = out[out["reactor"] == reactor]
            r.index = range(len(r.index))
            if od_filter:
                to_pop = []
                for i, OD, time in zip(r.index, r["od_measured"], r["exp_time"]):
                    if OD < od_filter:
                        to_pop.append(i)
                for i in to_pop:
                    r = r.drop(i)
                r.index = range(len(r.index))
            tmp = pd.DataFrame(columns=r.columns[:-1])
            j = 0
            while j + window_size < len(r):
                slice = r[j : j + window_size]
                row = [slice.iloc[-1]["exp_time"]]
                for c in slice.columns[1:-1]:
                    if c == "reactor":
                        row.append(reactor)
                    else:
                        row.append(np.average(slice[c]))
                tmp.loc[len(tmp)] = row
                j += window_size
            slice = r[j:]
            row = [slice.iloc[-1]["exp_time"]]
            for c in slice.columns[1:-1]:
                if c == "reactor":
                    row.append(reactor)
                else:
                    row.append(np.average(slice[c]))
            tmp.loc[len(tmp)] = row
            dfs.append(tmp)
        out = pd.concat(dfs)
    return out
