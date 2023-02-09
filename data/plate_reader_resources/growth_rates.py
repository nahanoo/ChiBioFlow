import numpy as np
import pandas as pd
import plotly.express as px
from math import log
df = pd.read_csv('data.csv')
df.insert(len(df.columns), 'dNdt', None)
df.index = range(len(df))
wells = set(df['well'])
for well in wells:
    tmp = df[df['well'] == well]
    dNdt = np.gradient(tmp['OD'], 0.5)/tmp['OD']
    for i,j in zip(tmp.index,dNdt):
        df.at[i,'dNdt'] = j
fig = px.line(df, x='Time', y='dNdt', line_group='well', color='species',
              line_dash='resource')
