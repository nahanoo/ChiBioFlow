import plotly.express as px
from model_species import Chain
import pandas as pd
from plotting import plot_species_model


chain = Chain(30)
chain.dilution_rate = 0.26
chain.transfer_rate = 1
chain.experiment(100)

"""for c in chain.chain[1:]:
    for specie in c.species.values():
        specie.N = 0
df = pd.DataFrame(columns=['x','N','species','reactor'])

i = 0
for c in chain.chain:
    for name,specie in c.species.items():
        for x,y in zip(specie.xs,specie.ys):
            df.loc[i] = [x,y,name,c.name]
            i += 1

fig = px.line(df,x='x',y='N',facet_col='reactor',color='species')
fig.show()"""
fig = plot_species_model(chain)
fig.show()