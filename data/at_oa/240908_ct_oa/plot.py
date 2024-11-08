import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, cfu_parser, ct_oa_parser
from plotly.subplots import make_subplots
from style_plot import *


def plot_exp():
    Ms, cfus = ct_oa_parser('/home/eric/ChiBioFlow/data/at_oa/240908_ct_oa')
    fig = make_subplots(rows=3,
                        cols=1,
                        shared_xaxes=True,
                        subplot_titles=[
                            'Total OD', 'Absolute abundances',
                            'Relative abundances'
                        ],
                        vertical_spacing=0.03)
    for M in Ms[:-1]:
        fig.add_trace(go.Scatter(x=M['exp_time'],
                                 y=M['od_measured'],
                                 marker=dict(color='#4F4F4F'),
                                 showlegend=False,
                                 hovertext=M['reactor']),
                      row=1,
                      col=1)
    for species in ['ct', 'oa']:
        for M in cfus[species][:-1]:
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
        for M in cfus[species][:-1]:
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
