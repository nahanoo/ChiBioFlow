from turtle import title
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go
from os.path import join


ms = ['m0','m1','m3','m4','m5','m8']

def mono_exp(x,m,t):
    return m * np.exp(-t * x)

def get_od(y,m,t):
    return np.log(y/m)/-t

def parse_df(f):
    return pd.read_csv(f,index_col=False)

def get_fit(f):
    df = parse_df(f)
    xs = 3*df['od'].to_list()
    ys = []
    ys += df['raw_1'].to_list()
    ys += df['raw_2'].to_list()
    ys += df['raw_3'].to_list()

    (m,t),cv = curve_fit(mono_exp,xs,ys)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xs,y=ys,mode='markers'))
    #fig.add_trace(go.Scatter(x=df['od'],y=mono_exp(df['od'],m,t),mode='lines'))
    fig.add_trace(go.Scatter(x=np.arange(0,5.1,0.05),y=mono_exp(np.arange(0,5.1,0.05),m,t),mode='lines'))
    fig.update_layout(
        title=f
    )
    fig.show()

for m in ms:
    f = join('calibrations',m+'.csv')
    get_fit(f)