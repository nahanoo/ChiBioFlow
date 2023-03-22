import pandas as pd
import plotly.express as px

cfus = pd.DataFrame(columns=['Treatment', 'species', 'CFUs'])
cfus.loc[0] = ['Mono D=0.1', 'O. anthropi', 0]
cfus.loc[1] = ['Mono D=0.1', 'C. Testosteroni', 2.53E9]
"""cfus.loc[2] = ['Coculture D=0', 'C. Testosteroni', 9.0E9]
cfus.loc[3] = ['Co D=0', 'O. anthropi', 9.3E9]"""
cfus.loc[2] = ['Co D=0.1', 'C. Testosteroni', 2.36E9]
cfus.loc[3] = ['Co D=0.1', 'O. anthropi', 3.1E9]
cfus.loc[4] = ['Co D=0.15', 'C. Testosteroni', 3.3E9]
cfus.loc[5] = ['Co D=0.15', 'O. anthropi', 7.43E8]


fig = px.scatter(cfus, x='Treatment', y='CFUs', color='species',
                 category_orders={'Treatment': ['Monoculture', 'Coculture']})
fig.update_xaxes(type='category')
fig.write_image('cfus.png', height=300, width=310)

citrate = pd.DataFrame(columns=['Treatment', 'Concentration [$\mu M$]'])
citrate.loc[0] = ['<i>C. testosteroni</i>', 2174]
citrate.loc[1] = ['Community', 1727]
fig = px.scatter(citrate, x='Treatment', y='Concentration [$\mu M$]', category_orders={
                 'Treatment': ['Monoculture', 'Coculture']})
fig.update_traces(marker=dict(color='#372f55'))
fig.update_xaxes(type='category')
fig.update_layout(title='Citric acid difference<br>feed steady state')
fig.write_image('conc.png', height=300, width=250)
