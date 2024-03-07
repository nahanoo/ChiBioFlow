import pandas as pd
import plotly.express as px
import numpy as np
from matplotlib import pyplot as plt

def dump_dfs():
    at = pd.read_csv('at.tsv',sep='\t')
    at.insert(len(at.columns),'name','at')
    oa = pd.read_csv('oa.tsv',sep='\t')
    oa.insert(len(oa.columns),'name','oa')

    df = pd.concat([at,oa])
    dfs = []

    for i,c in enumerate(df.columns[1:-1]):
        if i%3 == 0:
            carbon_source = i
        tmp = df[['Time',c,'name']]
        out = pd.DataFrame()
        out['time'],out['OD'],out['well'],out['strain'],out['carbon_source'] = tmp['Time'],tmp[c],c,tmp['name'],carbon_source
        dfs.append(out)

    df = pd.concat(dfs)
    df.to_csv('carbon_sources.csv',index=False)
    fum = df[df['carbon_source'] == 21]

    at = fum[fum['strain'] == 'at']
    at_rates = []
    for well in set(at['well']):
        OD = at[at['well'] == well]['OD']
        at_rates += list(np.gradient(OD) / OD * 2)
    at.insert(len(at.columns),'growth_rates',at_rates)

    oa = fum[fum['strain'] == 'oa']
    oa_rates = []
    for well in set(oa['well']):
        OD = oa[oa['well'] == well]['OD']
        oa_rates += list(np.gradient(OD) / OD * 2)
    oa.insert(len(oa.columns),'growth_rates',oa_rates)
    pd.concat([at,oa]).to_csv('fumarate.csv',index=False)


"""
fig = px.line(df,x='time',y='OD',facet_col='carbon_source',color='strain',facet_col_wrap=6,line_group='well')
fum = df[df['lg'] == 21]
fum_fig = px.line(fum,x='time',y='OD',color='strain',line_group='well')
"""