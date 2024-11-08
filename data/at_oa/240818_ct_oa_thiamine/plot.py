import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser
from plotly.subplots import make_subplots
from style_plot import *

# Data parsing
"""
dirs = [
    "240818_ct_oa_thiamine",
    "240818_ct_oa_thiamine_restart",
    "240818_ct_oa_thiamine_restart_2",
]

out_dir = "240818_ct_oa_thiamine_merged"

reactors = ["M0", "M1", "M2"]

for r in reactors:
    dfs = []
    for d in dirs:
        f = glob(
            join("/", "home", "eric", "ChiBioFlow", "data", d, r,
                 "*_data.csv"))[0]
        df = pd.read_csv(f)
        dfs.append(df)
    out = pd.concat(dfs)
    out["exp_time"] = np.linspace(0, 60 * len(out), len(out))
    out.to_csv(
        join("/", "home", "eric", "ChiBioFlow", "data", out_dir, r,
             split(f)[-1]))
"""


def plot_exp():
    Ms, cfus = ct_oa_parser(
        '/home/eric/ChiBioFlow/data/at_oa/240818_ct_oa_thiamine')

    fig = make_subplots(rows=3,
                        cols=1,
                        shared_xaxes=True,
                        subplot_titles=[
                            'Total OD', 'Absolute abundances',
                            'Relative abundances'
                        ],
                        vertical_spacing=0.03)
    for M in Ms:
        fig.add_trace(go.Scatter(x=M['exp_time'],
                                 y=M['od_measured'],
                                 marker=dict(color='#4F4F4F'),
                                 showlegend=False,
                                 hovertext=M['reactor']),
                      row=1,
                      col=1)
    for species in ['ct', 'oa']:
        for M in cfus[species]:
            fig.add_trace(go.Scatter(x=M['sample_time'],
                                     y=M['average'],
                                     name=abb[species],
                                     marker=dict(color=colors[species]),
                                     error_y=dict(type="data",
                                                  array=M["stdev"].to_list(),
                                                  visible=True),
                                     showlegend=False,
                                     hovertext=M['reactor']),
                          row=2,
                          col=1)

    for species in ['ct', 'oa']:
        for M in cfus[species]:
            fig.add_trace(go.Scatter(x=M['sample_time'],
                                     y=M['composition'],
                                     marker=dict(color=colors[species]),
                                     showlegend=False,
                                     hovertext=M['reactor']),
                          row=3,
                          col=1)

    fig['data'][3]['showlegend'] = True
    fig['data'][6]['showlegend'] = True
    fig["layout"]['yaxis2']["type"] = "log"
    fig["layout"]['yaxis1']["title"] = "OD"
    fig["layout"]['yaxis2']["title"] = "CFUs/mL"
    fig["layout"]['xaxis3']["title"] = "Time [h]"
    fig.update_layout(
        legend=dict(
            title='Species',
            y=0.6,  # Position the legend vertically; adjust as needed for placement
        ),
        width=width,
        height=height)
    fig = style_plot(fig)
    fig.show()
