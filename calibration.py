import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go
from os.path import join
import csv
import pandas as pd

reactors = ['M'+str(j) for j in range(8)]

def exp_fit(x,m,t):
    return m * np.exp(-t * x)

def dump_fit(plot=False):
    for reactor in reactors:
        fin = join('calibrations',reactor+'_raw.csv')
        xs = []
        ys = []
        df = pd.read_csv(fin,index_col=False)
        xs += df['od'].to_list() * 3
        for col in df.columns[1:]:
            ys += df[col].to_list()
        (m,t),cv = curve_fit(exp_fit,xs,ys)
        fout = join('calibrations',reactor+'_params.csv')
        df = pd.DataFrame(columns=['m','t'],index=[0])
        df.at[0,'m'] = m
        df.at[0,'t'] = t
        df.to_csv(fout,index=False)

        if plot:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=xs,y=ys,mode='markers'))
            #fig.add_trace(go.Scatter(x=df['od'],y=mono_exp(df['od'],m,t),mode='lines'))
            fig.add_trace(go.Scatter(x=np.arange(0,5.1,0.05),y=exp_fit(np.arange(0,5.1,0.05),m,t),mode='lines'))
            fig.update_layout(
                title=reactor
            )
            fig.show()
