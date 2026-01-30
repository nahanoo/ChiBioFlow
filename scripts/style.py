width, height = 150, 150
colors = {
    "blue": "#000080",
    "ct": "#7570B3",
    "oa": "#D95F02",
    "ms": "#E6AB02",
    "at": "#1B9E77",
    "Spent media Ct": "#1B9E77",
    "Spent media Oa": "#E7298A",
    "H20": "gray",
}

colors_heatmap = [
    [0.0, "#7570B3"],  # deep blue
    [0.5, "white"],
    [1.0, "#D95F02"],  # deep red
]

colors_metabolites = {
    "Nucleotide related": "#1f77b4",
    "Carbohydrates": "#ff7f0e",
    "Fatty Acids": "#2ca02c",
    "Amino Acids": "#9467bd",
    "Organic Acids": "#8c564b",
    "Coenzymes": "#e377c2",
    "Others": "#7f7f7f",
}


def style_plot(
    fig,
    marker_size=3,
    top_margin=10,
    left_margin=10,
    right_margin=10,
    buttom_margin=10,
    font_size=14,
    line_thickness=1.5,
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
    gridline_width = 0.2
    fig.update_yaxes(
        title_standoff=0,
        gridcolor="gray",
        zeroline=False,
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=0.5,
        showline=True,
        mirror=True,
        linecolor="black",
        linewidth=0.5,
        tickcolor="black",
        tickwidth=0.5,
    )
    fig.update_xaxes(
        title_standoff=0,
        gridcolor="gray",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=0.5,
        showline=True,
        mirror=True,
        linecolor="black",
        linewidth=0.5,
        zeroline=False,
        tickcolor="black",
        tickwidth=0.5,
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    return fig


def square_panel_by_height(
    fig,
    height_px: int,
    cbar_trace: int = 0,
    cbar_gap_px: int = 8,
    cbar_label_px: int = 40,  # budget for tick labels (and title if you add one)
    extra_right_px: int = 0,  # any extra whitespace you want
):
    # 1) total figure height is fixed
    fig.update_layout(autosize=False, height=height_px)

    # margins (must be set before calling this; e.g. by your style_plot)
    m = fig.layout.margin
    l = int(m.l or 0)
    r = int(m.r or 0)
    t = int(m.t or 0)
    b = int(m.b or 0)

    # 2) square panel side = inner height
    side = height_px - t - b
    side = max(side, 10)  # safety

    # 3) reserve right-side pixels for colorbar + labels
    thick = 0
    if cbar_trace is not None:
        cb = getattr(fig.data[cbar_trace], "colorbar", None)
        thick = int(getattr(cb, "thickness", 0) or 0)

    right_px = thick + cbar_gap_px + cbar_label_px + extra_right_px

    inner_w = side + right_px
    fig.update_layout(width=l + r + inner_w)

    # 4) shrink x-domain so plot area is exactly `side` pixels wide
    dom_end = side / inner_w
    fig.update_xaxes(domain=[0.0, dom_end])
    fig.update_yaxes(domain=[0.0, 1.0])

    # 5) put the colorbar into the reserved right-side area (still in "paper" coords) :contentReference[oaicite:1]{index=1}
    if cbar_trace is not None and thick > 0:
        x_left = dom_end + (cbar_gap_px / inner_w)
        fig.data[cbar_trace].update(colorbar=dict(x=x_left, xanchor="left"))

    return fig
