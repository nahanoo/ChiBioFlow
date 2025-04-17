colors = {
    "ct": "#7570B3",
    "oa": "#D95F02",
    "<i>C. testosteroni</i>": "#7570B3",
    "<i>O. anthropi</i>": "#D95F02",
    "total": "#4F4F4F",
    "cartoon": "#5b4b8a",
}

abb = {
    "ct": "<i>C. testosteronis</i>",
    "oa": "<i>O. anthropi</i>",
    "<i>C. testosteroni</i>": "at",
    "<i>O. anthropi</i>": "oa",
}

width, height = 1000, 800


def style_plot(
    fig,
    marker_size=3,
    top_margin=20,
    left_margin=40,
    right_margin=0,
    buttom_margin=50,
    font_size=14,
):
    """Style function for figures setting fot size and true black color."""
    fig.update_layout(
        {
            "plot_bgcolor": "white",
            "paper_bgcolor": "white",
        },
        font={"size": font_size, "color": "black"},
    )
    for d in fig["data"]:
        d["marker"]["size"] = marker_size
        d["line"]["width"] = marker_size
        d["error_y"]["thickness"] = marker_size
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = font_size
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = font_size
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = font_size
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"

    fig.update_layout(
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=buttom_margin),
        hoverlabel=dict(font_size=font_size),
    )
    gridline_width = 1
    fig.update_yaxes(
        title_standoff=0,
        gridcolor="black",
        # zeroline=True,
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
        # minor=dict(showgrid=False),
    )
    fig.update_xaxes(
        title_standoff=0,
        # zeroline=True,
        gridcolor="black",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
        # minor=dict(showgrid=False)
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    return fig
