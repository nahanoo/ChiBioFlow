import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser
from plotly.subplots import make_subplots
from style_plot import *


def plot_exp():
    Ms, cfus = ct_oa_parser('/home/eric/ChiBioFlow/data/at_oa/240921_ct',
                            parse_cfus=False)
    fig = make_subplots(rows=1,
                        cols=1,
                        shared_xaxes=True,
                        subplot_titles=['OD'])
    for M in Ms:
        fig.add_trace(go.Scatter(x=M['exp_time'],
                                 y=M['od_measured'],
                                 marker=dict(color='#4F4F4F'),
                                 showlegend=False,
                                 hovertext=M['reactor']),
                      row=1,
                      col=1)

    fig["layout"]['yaxis1']["title"] = "OD"
    fig["layout"]['xaxis1']["title"] = "Time [h]"
    fig.update_layout(width=width, height=height / 2)
    fig = style_plot(fig)
    fig.show()


plot_exp()
