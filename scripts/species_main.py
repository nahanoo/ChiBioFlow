from model_species import Chain
import pandas as pd
import numpy as np
import plotly.express as px


def random_compositions(t,n):
    dfs = []
    compositions = [np.random.dirichlet(np.ones(4), size=1)[
        0] for i in range(n)]
    for counter, composition in enumerate(compositions):
        chain = Chain(4)
        chain.dilution_rate = 0.2
        chain.transfer_rate = 2
        for j, specie in enumerate(chain.chain[0].species.values()):
            specie.N = np.random.choice(np.arange(0.08,2,0.2)) * compositions[counter][j]

        for c in chain.chain[1:]:
            for specie in c.species.values():
                specie.N = 0

        chain.experiment(t)

        for c in chain.chain:
            for name, specie in c.species.items():
                df = pd.DataFrame({'x': chain.xs, 'N': specie.ys,
                                   'species': name, 'reactor': c.name, 'composition': counter})
                dfs.append(df)
            

    out = pd.concat(dfs)
    fig = px.line(out, x='x', y='N', facet_col='reactor',
                color='species', line_group='composition')
    return fig

def random_dilution_rates():
    dfs = []
    dilution_fs = np.arange(0.05, 0.45, 0.05)
    for counter, dilution_f in enumerate(dilution_fs):
        chain = Chain(4)
        chain.dilution_rate = dilution_f
        chain.transfer_rate = 30
        for j, specie in enumerate(chain.chain[0].species.values()):
            specie.N = 0.02

        for c in chain.chain[1:]:
            for specie in c.species.values():
                specie.N = 0

        chain.experiment(100)

        for c in chain.chain:
            for name, specie in c.species.items():
                df = pd.DataFrame({'x': chain.xs, 'N': specie.ys,
                                   'species': name, 'reactor': c.name, 'dilution factor': dilution_f})
                dfs.append(df)
            df = pd.DataFrame({'x': chain.xs, 'N': c.total,
                               'species': 'total', 'reactor': c.name, 'dilution factor': dilution_f})
            dfs.append(df)

    out = pd.concat(dfs)
    f = px.line(out, x='x', y='N', facet_col='reactor',
                color='species', line_group='dilution factor')
    f.show()
