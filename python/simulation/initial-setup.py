#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
from argparse import ArgumentParser, Namespace
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_arguments() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Visualize simulation results.")
    parser.add_argument("--out-fig-dir", type=str, help="Output figure directory")

    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def draw_anomaly(
    model_type,
    ax,
    X,
    Y,
    T_anomaly,
    cmap_modified,
    n_colors,
    x_extent,
    y_extent,
    phase_transition_depth,
    text_offset,
    border_offset,
    arrow_offset,
    n_arrows,
    x_velocity,
    y_velocity,
    velocity_factor,
    arrow_head_width,
    arrow_head_length,
    bc_labels: dict,
    show_axes: bool,
):
    """"""
    levels = np.linspace(-500, 500, n_colors + 1)
    contour = ax.contourf(X, Y, T_anomaly, cmap=cmap_modified, levels=levels, vmin=-500, vmax=500)

    # phase transition
    ax.hlines(y=phase_transition_depth, xmin=0, xmax=x_extent, color="black", linewidth=0.5)
    ax.text(x_extent / 2, phase_transition_depth - text_offset, "Olivine", va="bottom", ha="center")
    ax.text(x_extent / 2, phase_transition_depth + text_offset, "Wadsleyite", va="top", ha="center")

    # velocity boundary conditions (arrows)
    x_positions_top = np.linspace(arrow_offset * 2, x_extent - (arrow_offset * 2), n_arrows)
    y_positions_right = np.linspace(0 + (3 * arrow_offset), y_extent - arrow_offset, n_arrows)
    x_positions_bottom = np.linspace(arrow_offset * 2, x_extent - (arrow_offset * 2), n_arrows)
    y_positions_left = np.linspace(0 + (3 * arrow_offset), y_extent - arrow_offset, n_arrows)

    if model_type == "slab":
        for x in x_positions_top:
            ax.arrow(
                x,
                0,
                x_velocity * velocity_factor,
                y_velocity * velocity_factor,
                head_width=arrow_head_width,
                head_length=arrow_head_length,
                color="black",
            )
    if model_type == "slab":
        for y in y_positions_right:
            ax.arrow(x_extent - 10, y, x_velocity * velocity_factor, 0, head_width=arrow_head_width, head_length=arrow_head_length, color="black")
    if model_type == "slab":
        for x in x_positions_bottom:
            ax.arrow(
                x,
                y_extent,
                x_velocity * velocity_factor,
                0,
                head_width=arrow_head_width,
                head_length=arrow_head_length,
                color="black",
            )
    else:
        for x in x_positions_bottom:
            ax.arrow(
                x,
                y_extent,
                x_velocity * velocity_factor,
                y_velocity * velocity_factor,
                head_width=arrow_head_width,
                head_length=arrow_head_length,
                color="black",
            )
    if model_type == "slab":
        for y in y_positions_left:
            ax.arrow(0, y, x_velocity * velocity_factor, 0, head_width=arrow_head_width, head_length=arrow_head_length, color="black")

    # boundary condition labels
    ax.text(x_extent / 2, 0 - text_offset, bc_labels["top"], va="bottom", ha="center")
    ax.text(x_extent + text_offset, y_extent / 2, bc_labels["right"], va="center", ha="left", rotation=270)
    ax.text(-text_offset, y_extent / 2, bc_labels["left"], va="center", ha="right", rotation=90)
    ax.text(x_extent / 2, y_extent + text_offset, bc_labels["bottom"], va="top", ha="center")

    # formatting
    ax.set_xlim(-border_offset, x_extent + border_offset)
    ax.set_ylim(-border_offset, y_extent + border_offset)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Depth (km)")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if not show_axes:
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

    return contour


#######################################################
## .1. Main                                      !!! ##
#######################################################
def main():
    """"""
    args = parse_arguments()

    out_fig_dir = Path(args.out_fig_dir) if args.out_fig_dir else Path("./figures")
    if not out_fig_dir.exists():
        out_fig_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_fig_dir / "initial-setup.png"

    if out_path.exists():
        print(f" -- Found plot: {out_path.name}!")
        return

    print(f" -> {out_path.name}")

    plot_width: float = 5.0
    plot_height: float = 7.0
    n_colors: float = 13
    central_color_hex: str = "#E6E6E6"

    x_extent: float = 396
    y_extent: float = 264
    phase_transition_depth: float = 132

    anomaly_center_slab: tuple[float, float] = (118, 0)
    anomaly_sigma_slab: float = 15
    anomaly_amplitude_slab: float = -500
    x_velocity_slab: float = 0.5
    y_velocity_slab: float = 1.5

    anomaly_center_plume: tuple[float, float] = (x_extent / 2, y_extent)
    anomaly_sigma_plume: float = 15
    anomaly_amplitude_plume: float = 500
    x_velocity_plume: float = 0
    y_velocity_plume: float = -1.5

    text_offset: float = 10
    border_offset: float = 40
    n_arrows: float = 7
    arrow_offset: float = 10
    arrow_head_width: float = 10
    arrow_head_length: float = 5
    velocity_factor: float = 10

    rcParams = {
        "figure.dpi": 300,
        "savefig.bbox": "tight",
        "axes.facecolor": "1.0",
        "legend.frameon": False,
        "legend.facecolor": "0.9",
        "legend.loc": "upper left",
        "legend.fontsize": "small",
        "figure.autolayout": True,
        "font.size": 12,
    }

    nx, ny = 300, 200
    x, y = np.linspace(0, x_extent, nx), np.linspace(y_extent, 0, ny)
    X, Y = np.meshgrid(x, y)
    x0_slab, y0_slab = anomaly_center_slab
    x0_plume, y0_plume = anomaly_center_plume
    T_anomaly_slab = anomaly_amplitude_slab * np.exp(-(((Y - y0_slab) ** 2 + (X - x0_slab) ** 2) / (2 * anomaly_sigma_slab**2)))
    T_anomaly_plume = anomaly_amplitude_plume * np.exp(-(((Y - y0_plume) ** 2 + (X - x0_plume) ** 2) / (2 * anomaly_sigma_plume**2)))

    original_cmap = plt.get_cmap("seismic")
    cmap_colors = original_cmap(np.linspace(0, 1, n_colors))
    center_idx = n_colors // 2
    cmap_colors[center_idx] = mcolors.to_rgba(central_color_hex)
    cmap_modified = mcolors.ListedColormap(cmap_colors)

    fig = plt.figure(figsize=(plot_width, plot_height), constrained_layout=True)
    plt.rcParams.update(rcParams)
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.05])

    ax_slab = fig.add_subplot(gs[0, 0])
    ax_plume = fig.add_subplot(gs[1, 0])

    # Define different BC labels for slabs vs plumes
    bc_slab = {
        "top": "Prescribed inflow",
        "right": "$\\sigma_{xy}$ = 0",
        "left": "$\\sigma_{xy}$ = 0",
        "bottom": "$\\sigma_{yy}$ = $\\bar{P}$",
    }

    bc_plume = {
        "top": "$\\sigma_{yy}$ = $\\bar{P}$",
        "right": "$\\sigma_{xy}$ = 0",
        "left": "$\\sigma_{xy}$ = 0",
        "bottom": "Prescribed inflow",
    }

    # Draw slab anomaly
    contour_slab = draw_anomaly(
        "slab",
        ax_slab,
        X,
        Y,
        T_anomaly_slab,
        cmap_modified,
        n_colors,
        x_extent,
        y_extent,
        phase_transition_depth,
        text_offset,
        border_offset,
        arrow_offset,
        n_arrows,
        x_velocity_slab,
        y_velocity_slab,
        velocity_factor,
        arrow_head_width,
        arrow_head_length,
        bc_slab,
        False,
    )

    # Draw plume anomaly
    _ = draw_anomaly(
        "plume",
        ax_plume,
        X,
        Y,
        T_anomaly_plume,
        cmap_modified,
        n_colors,
        x_extent,
        y_extent,
        phase_transition_depth,
        text_offset,
        border_offset,
        arrow_offset,
        n_arrows,
        x_velocity_plume,
        y_velocity_plume,
        velocity_factor,
        arrow_head_width,
        arrow_head_length,
        bc_plume,
        True,
    )

    # Shared colorbar
    ticks = [-500, -250, 0, 250, 500]
    cax = fig.add_axes([0.34, 0.43, 0.4, 0.015])  # pyright: ignore
    cbar = fig.colorbar(contour_slab, cax=cax, ticks=ticks, orientation="horizontal")
    cbar.set_label("Thermal Anomaly (K)", fontsize="11")
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.xaxis.set_label_position("top")

    plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.05)


if __name__ == "__main__":
    main()
