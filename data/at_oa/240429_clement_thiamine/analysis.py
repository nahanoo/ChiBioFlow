import pandas as pd
import plotly.express as px


def parse_data():
    dfs = []
    data = pd.read_csv('thiamine.txt')
    time = data['Time']
    for c in data.columns[2:]:
        df = pd.DataFrame(columns=['time','OD','well'])
        df['time'],df['OD'],df['well'] = time, data[c], c
        dfs.append(df)

    out = pd.concat(dfs)
    out.to_csv('plate_reader.csv',index=False)

df = pd.read_csv('plate_reader.csv')

def plot_wells():
    fig = px.line(df,x='time',y='OD',facet_col='well',facet_col_wrap=12)
    fig.show()
    fig.write_html('thiamine.html')

thiamine_stock = 0.025 / 337.27 / 0.1 / 500 * 1000
B4_t0 = thiamine_stock * 0.75 * (0.25 **3)
B4 = df[df['well'] == 'B4']
B4_N1 = 0.51
B4_Y = B4_N1 / B4_t0

C4 = df[df['well'] == 'C4']
C4_t0 = thiamine_stock * 0.5 * (0.25 **3)
C4_N1 = 0.41
C4_Y = C4_N1 / C4_t0
#px.line(B4,x='time',y='OD',facet_col='well',facet_col_wrap=12).show()
#px.line(C4,x='time',y='OD',facet_col='well',facet_col_wrap=12).show()

