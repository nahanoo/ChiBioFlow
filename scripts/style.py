width, height = 150, 150
colors = {
    "blue": "#000080",
    "ct": "#7570B3",
    "oa": "#D95F02",
    "Spent media Ct": "#1B9E77",
    "Spent media Oa": "#E7298A",
    "H20": "gray",
}


def style_plot(
    fig,
    marker_size=3,
    top_margin=10,
    left_margin=10,
    right_margin=10,
    buttom_margin=10,
    font_size=14,
    line_thickness=3,
):
    """Style function for figures setting fot size and true black color."""
    fig.update_layout(
        {
            "plot_bgcolor": "#FFFFFF",
            "paper_bgcolor": "#FFFFFF",
        },
        font={"size": font_size, "color": "black"},
    )
    for d in fig["data"]:
        try:
            d["marker"]["size"] = marker_size
        except KeyError:
            pass
        try:
            d["line"]["width"] = line_thickness
        except KeyError:
            pass
        try:
            d["error_y"]["thickness"] = line_thickness
        except KeyError:
            pass
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
    gridline_width = 0.5
    fig.update_yaxes(
        title_standoff=0,
        gridcolor="black",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
    )
    fig.update_xaxes(
        title_standoff=0,
        gridcolor="black",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    return fig
