import numpy as np
import pandas as pd
import plotly.express as px

df = pd.read_csv('data.csv')
df.insert(len(df.columns), 'dNdt', None)
df['dNdt'] = np.gradient(df['OD'], 0.5)/df['OD']
df = df[df['Time'] >= 1]
df = df[df['Time'] <= 74]
fig = px.line(df, x='Time', y='dNdt', line_group='well', color='species',
              line_dash='resource')
