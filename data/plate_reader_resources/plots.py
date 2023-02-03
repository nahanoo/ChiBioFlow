import plotly.express as px
import pandas as pd
import numpy as np

df = pd.read_csv('ODs_export.csv')
cs = pd.read_csv('conditions.csv')
cs.index = cs['well']

dfs = []
columns = cs.index
for c in columns:
    out = pd.DataFrame(columns=['Time', 'OD', 'well', 'resource', 'species'])
    out['Time'] = np.arange(0.5,76,0.5)
    out['OD'] = df[c]
    out['well'] = c
    out['resource'] = cs.loc[c]['resource']
    out['species'] = cs.loc[c]['species']
    dfs.append(out)
out = pd.concat(dfs)
out.to_csv('data.csv',index=False)
fig = px.line(out, x='Time', y='OD', line_group='well', color='species',
              line_dash='resource')
