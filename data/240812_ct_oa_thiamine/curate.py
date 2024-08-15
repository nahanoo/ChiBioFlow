import pandas as pd
from os.path import join, split
import numpy as np
from glob import glob
from subprocess import call
from os import symlink
import plotly.express as px

reactors = ["M0", "M1", "M2"]


def mask_data(t1, t2):
    for r in reactors:
        f_in = glob(join(r, "*data.csv"))[0]
        df = pd.read_csv(f_in)
        for i, row in df.iterrows():
            if row["exp_time"] >= t1 * 60**2:
                t0 = row["exp_time"]
                break
        df = df.loc[i:]
        df.loc[:, "exp_time"] = df["exp_time"] - t0
        for i, row in df.iterrows():
            if row["exp_time"] >= t2 * 60**2:
                break
        df = df.loc[:i]
        f_out = join(
            "../240812_ct_oa_thiamine_curated",
            r,
            split(f_in)[1],
        )
        df.to_csv(f_out, index=False)


def plot():
    cmd = ["python", "~/ChiBioFlow/data/at_oa/240812_ct_oa_thiamine/plot.py"]
    call(" ".join(cmd), shell=True)


def fix_M1():
    f = glob(join("../240812_ct_oa_thiamine_curated", "M1", "*data.csv"))[0]
    df = pd.read_csv(f)
    df.loc[541:551, "od_measured"] = 0.416
    df.loc[:540, "od_measured"] += 0.215
    df = df.drop(2350)
    df.loc[:, "od_measured"] = df.loc[:, "od_measured"] - 0.0683
    df.loc[:538, "FP1_base"] += 15000
    df.loc[:538, "FP1_emit1"] -= 0.0666453674121405 - 0.0564727507914354
    i = 20
    df.loc[538 : 538 + i, "FP1_base"] = df.loc[538, "FP1_base"]
    df.loc[538 : 538 + i, "FP1_emit1"] = df.loc[538, "FP1_emit1"]
    df.to_csv(f, index=False)


def fix_M1():
    f = glob(join("../240812_ct_oa_thiamine_curated", "M0", "*data.csv"))[0]
    df = pd.read_csv(f)
    for i, row in df[df["od_measured"] < 0.15].iterrows():
        df = df.drop(i)
    df.loc[:, "od_measured"] -= 0.03
    # px.line(df, x="exp_time", y="od_measured").show()
    df.to_csv(f)


def fix_M2():
    f = glob(join("../240812_ct_oa_thiamine_curated", "M2", "*data.csv"))[0]
    df = pd.read_csv(f)
    df.loc[1388:1402, "od_measured"] = 0.275
    df.loc[:1387, "od_measured"] += 0.172
    df = df.drop(2351)
    df.loc[2352:, "od_measured"] += 0.005
    df.loc[:1388, "FP1_base"] += 13000
    df.loc[:1388, "FP1_emit1"] -= 0.113002862739189 - 0.0944862259525491 + 0.004
    i = 20
    df.loc[1388 : 1388 + i, "FP1_base"] = df.loc[1388, "FP1_base"]
    df.loc[1388 : 1388 + i, "FP1_emit1"] = df.loc[1388, "FP1_emit1"]
    df.to_csv(f)


mask_data(1.71, 63.08)
fix_M1()
fix_M2()

plot()

"""f = glob(join("../240812_ct_oa_thiamine_curated", "M1", "*data.csv"))[0]
df = pd.read_csv(f)
df.loc[:541, "od_measured"] = df.loc[:541, "od_measured"] + 0.215
df.loc[:, "od_measured"] = df.loc[:, "od_measured"] - 0.0683
f_out = join(
    "~/ChiBioFlow/data/240812_ct_oa_thiamine_curated",
    "M1",
    split(f)[1].replace(".xlsx", ".csv"),
)
df.to_csv(f_out, index=False)
"""
