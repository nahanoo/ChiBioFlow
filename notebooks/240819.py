import plotly.graph_objects as go
import numpy as np
from glob import glob
import pandas as pd
from fitting import *


meta, raw = get_dfs("./")
# Filter data
mask = (meta["species"] == "Comamonas testosteroni") & (meta["cs_conc"] == 15)
columns = list(set(meta[mask]["linegroup"]))
ct15 = raw[["time"] + columns].dropna()
