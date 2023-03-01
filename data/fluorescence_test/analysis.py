import pandas as pd
from glob import glob
import plotly.express as px
import numpy as np

wells = ['Time', 'B2', 'B3', 'B4', 'C2', 'C3', 'C4',
         'D2', 'D3', 'D4', 'E2', 'E3', 'E4']
dfs = []
fs = glob('*_*.csv')
for f in fs:
    name = f.split('_')[0]
    df = pd.read_csv(f, usecols=wells)
    df.insert(0, 'fluorophore', name)
    dfs.append(df)
    time = df['Time']
df = pd.concat(dfs)

meta = pd.read_csv('conditions.csv')
meta.index = meta['well']
dfs = []
columns = meta.index
fluorophores = ['gfp', 'morange', 'mcherry']
for fluorophore in fluorophores:
    tmp = df[df['fluorophore'] == fluorophore]
    for c in columns:
        out = pd.DataFrame(
            columns=['Time', 'intensity', 'well', 'species', 'fluorophore'])
        out['Time'] = time
        out['intensity'] = tmp[c]
        out['well'] = c
        out['species'] = meta.loc[c]['species']
        out['fluorophore'] = tmp['fluorophore']
        print(out)
        dfs.append(out)
data = pd.concat(dfs)
fig = px.line(data, x='Time', y='intensity', line_group='well',
              color='species',line_dash='fluorophore')
fig.write_html('fluorescence.html')