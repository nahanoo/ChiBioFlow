import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from model_class import Chemostat
import plotly.express as px

c1 = Chemostat('data/overnight_08_05/M6/2021-04-26 21_45_59_M6_data.csv',0.05)
chain = [c1]

c1.plot_fit()
dilution = c1.get_dilution()
print(dilution)
def experiment():
    x0 = 0
    for i in range(48):
        for counter,c in enumerate(chain):
            t0 = c.get_x(c.N0)
            t1 = t0 + 3600
            c.N1 = c.Nt(t1,c.N0)
            if counter == 0:
                c.N0 = c.N1 * 20 / (20 + dilution)
            else:
                c.N0 = (chain[counter - 1].N1 * dilution + c.N1 * 20) / (dilution + 20)
            xs = np.linspace(x0,x0 + 3600, 6)
            ys = [c.Nt(t,c.N0) for t in np.linspace(t0,t1,6)]
            c.xs = np.concatenate([c.xs,xs])
            c.ys = np.concatenate([c.ys,ys])
        x0 = x0 + 3600 +60

    df = pd.DataFrame(columns=['Time','OD','c'])
    for counter,c in enumerate(chain):
        tmp = pd.DataFrame(columns=['Time','OD','c'])
        tmp['Time'] = [float(e)/3600 for e in c.xs]
        tmp['OD'] = [float(e) for e in c.ys]
        tmp['c'] = len(c.xs) * [counter]
        df = pd.concat([df,tmp])

    f = px.line(df,x='Time',y='OD',facet_col='c',facet_col_wrap=2)
    f.show()