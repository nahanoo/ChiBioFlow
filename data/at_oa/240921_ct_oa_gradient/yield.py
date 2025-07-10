import pandas as pd
import numpy as np
from chibio_parser import *
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from os.path import join

species_dict = {"Comamonas testosteroni": "ct", "Ochrobactrum anthropi": "oa"}

dir = join("/", "home", "eric", "ChiBioFlow", "data")
cfus = cfu_parser(
    join(dir, "at_oa/240921_ct_oa_gradient"),
)[0]
df = pd.DataFrame(
    columns=["species", "concentration", "OD", "OD_error", "CFUs", "CFUs_error"]
)

meta_full = pd.read_csv("/home/eric/curves/metadata/pooled_df_joint_metadata.csv")
raw = pd.read_csv(
    "/home/eric/curves/export/240921_carbon_source_gradients/measurement_data.csv"
)
raw.index = raw["linegroup"]
for species in ["Comamonas testosteroni", "Ochrobactrum anthropi"]:
    mask = (meta_full["project"] == "ct_oa_chemostat_project/") & (
        meta_full["species"] == species
    )
    meta = meta_full[mask]
    for i, c in enumerate(set(meta["cs_conc"])):
        yield_ods = []
        c_sub = meta[meta["cs_conc"] == c]
        for l in c_sub["linegroup"]:
            x, y = (
                raw.loc[l, "time"].to_numpy(),
                raw.loc[l, "measurement"].to_numpy(),
            )
            yield_ods.append(max(y))
        cfu_sub = cfus[
            (cfus["species"] == species_dict[species]) & (cfus["comment"] == c)
        ]
        cfu_yield, cfu_yield_std = cfu_sub[["average", "stdev"]].iloc[0].to_list()

        df.loc[species.replace(" ", "_") + "_" + str(c)] = [
            species,
            c,
            np.average(np.array(yield_ods)),
            np.std(np.array(yield_ods)),
            cfu_yield,
            cfu_yield_std,
        ]

fig = make_subplots(
    rows=2, cols=2, row_titles=["Comamonas testosteroni", "Ochrobactrum anthropi"]
)
df = df.sort_values("concentration")
for i, species in enumerate(["Comamonas testosteroni", "Ochrobactrum anthropi"]):
    df_sub = df[df["species"] == species]
    fig.add_trace(
        go.Scatter(
            x=df_sub["concentration"],
            y=df_sub["OD"],
            error_y=dict(type="data", array=df_sub["OD_error"]),
        ),
        row=1 + i,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=df_sub["concentration"],
            y=df_sub["CFUs"],
            error_y=dict(type="data", array=df_sub["CFUs_error"]),
        ),
        row=1 + i,
        col=2,
    )
fig.update_yaxes(type="log", row=[1, 2], col=1)
fig.show()
