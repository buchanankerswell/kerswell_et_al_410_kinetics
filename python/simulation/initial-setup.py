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
    """"""
    parser = ArgumentParser(description="Visualize simulation results.")
    parser.add_argument("--out-fig-dir", type=str, help="Output figure directory")

    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_thermal_anomaly(X, Y, x_start, y_start, x_end, y_end, w, dT, alpha):
    """"""
    x_start = x_start
    y_start = y_start
    x_end = x_end
    y_end = y_end
    w = w
    alpha = alpha

    dx = x_end - x_start
    dy = y_end - y_start
    L = np.sqrt(dx**2 + dy**2)

    x_shifted = X - x_start
    y_shifted = Y - y_start

    perp_dist = x_shifted * (-dy) + y_shifted * dx
    parallel_dist = x_shifted * dx + y_shifted * dy

    gaussian_term = np.exp(-(perp_dist**2) / (2 * w**2 * L**2))
    tanh1 = 1 + np.tanh((parallel_dist / L) / alpha)
    tanh2 = 1 + np.tanh((L - parallel_dist / L) / alpha)

    T = dT * gaussian_term * 0.25 * tanh1 * tanh2

    return T


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def draw_initial_conditions(ax, X, Y, T, cmap, x_extent, y_extent, phase_depth, bc_labels, title, text_offset=25):
    """"""
    levels = np.linspace(-500, 500, 14)
    contour = ax.contourf(X, Y, T, levels=levels, cmap=cmap)

    ax.axhline(y=phase_depth, color="black", linestyle="--", linewidth=1.5, zorder=5)

    ax.set_xlim(0, x_extent)
    ax.set_ylim(0, y_extent)
    ax.set_aspect("equal")

    ax.text(
        x_extent - text_offset if "Slab" in title else x_extent / 2,
        y_extent - text_offset,
        bc_labels["top1"],
        ha="right" if "Slab" in title else "center",
        va="top",
        fontsize=13,
    )

    ax.text(
        text_offset if "Slab" in title else x_extent / 2,
        y_extent - text_offset,
        bc_labels["top2"],
        ha="left" if "Slab" in title else "center",
        va="top",
        fontsize=13,
    )

    ax.text(
        x_extent / 2 if "Slab" in title else text_offset,
        text_offset,
        bc_labels["bottom1"],
        ha="center" if "Slab" in title else "left",
        va="bottom",
        fontsize=13,
    )

    ax.text(
        x_extent / 2 if "Slab" in title else x_extent - text_offset,
        text_offset,
        bc_labels["bottom2"],
        ha="center" if "Slab" in title else "right",
        va="bottom",
        fontsize=13,
    )

    ax.text(text_offset * 2.6, y_extent / 2, bc_labels["left"], ha="center", va="center", rotation=90, fontsize=13)
    ax.text(x_extent - text_offset * 2.6, y_extent / 2, bc_labels["right"], ha="center", va="center", rotation=90, fontsize=13)

    if "Slab" in title:
        ax.set_xticks([])

    if "Slab" not in title:
        ax.set_xlabel("X (km)")

    ax.set_ylabel("Y (km)")
    ax.set_title(title, fontsize=16)

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

    print("    --------------------------------------------------")
    print("    Drawing initial setup")
    print("    --------------------------------------------------")

    if out_path.exists():
        print(f" -- Found plot: {out_path.name}!")
        return

    print(f" -> {out_path.name}")

    plot_width = 4.5
    plot_height = 6.5
    n_colors = 13
    central_color_hex = "#E6E6E6"

    x_extent = 900
    y_extent = 600
    phase_transition_depth_slab = 473
    phase_transition_depth_plume = 460

    nx, ny = 900, 600
    x = np.linspace(0, x_extent, nx)
    y = np.linspace(0, y_extent, ny)
    X, Y = np.meshgrid(x, y)

    slab_params = {
        "x_start": 350,
        "y_start": 630,
        "x_end": 450,
        "y_end": 500,
        "w": 15,
        "dT": -500,
        "alpha": 5,
        "dx": 100,
        "dy": -130,
        "L": np.sqrt(100**2 + 130**2),
        "v": 0.05,
    }

    plume_params = {
        "x_start": 450,
        "y_start": -30,
        "x_end": 450,
        "y_end": 450,
        "w": 15,
        "dT": 500,
        "alpha": 5,
        "dx": 0,
        "dy": 450,
        "L": np.sqrt(0**2 + 450**2),
        "v": 0.05,
    }

    T_slab = compute_thermal_anomaly(X, Y, **{k: slab_params[k] for k in ["x_start", "y_start", "x_end", "y_end", "w", "dT", "alpha"]})
    T_plume = compute_thermal_anomaly(X, Y, **{k: plume_params[k] for k in ["x_start", "y_start", "x_end", "y_end", "w", "dT", "alpha"]})

    n_colors = 13
    central_color_hex = "#E6E6E6"
    original_cmap = plt.get_cmap("seismic")
    cmap_colors = original_cmap(np.linspace(0, 1, n_colors))
    center_idx = n_colors // 2
    cmap_colors[center_idx] = mcolors.to_rgba(central_color_hex)
    cmap_modified = mcolors.ListedColormap(cmap_colors)

    bc_slab = {
        "top1": "$\\vec{u}_x$=$f(x)$, $\\vec{u}_y$=$f(x)$",
        "top2": "Fixed $T$",
        "right": "$\\sigma_{xx}$ = $d\\bar{P}/dy$,\n $\\vec{u}_x$=0",
        "left": "$\\sigma_{xx}$ = $d\\bar{P}/dy$,\n $\\vec{u}_x$=0",
        "bottom1": "$\\sigma_{yy}$ = $\\bar{P}(bottom)$, $\\vec{u}_x$=0",
        "bottom2": None,
    }

    bc_plume = {
        "top1": "$\\sigma_{yy}$ = $\\bar{P}(top)$, $\\vec{u}_x$=0",
        "top2": None,
        "right": "$\\sigma_{xx}$ = $d\\bar{P}/dy$,\n $\\vec{u}_x$=0",
        "left": "$\\sigma_{xx}$ = $d\\bar{P}/dy$,\n $\\vec{u}_x$=0",
        "bottom1": "Fixed $T$",
        "bottom2": "$\\vec{u}_x$=0, $\\vec{u}_y$=$f(x)$",
    }

    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.bbox": "tight",
            "axes.facecolor": "0.9",
            "legend.frameon": False,
            "legend.facecolor": "0.9",
            "legend.loc": "upper left",
            "legend.fontsize": "small",
            "figure.autolayout": True,
            "font.size": 14,
        }
    )

    fig = plt.figure(figsize=(plot_width, plot_height), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.05])

    ax_slab = fig.add_subplot(gs[0, 0])
    ax_plume = fig.add_subplot(gs[1, 0])

    contour_slab = draw_initial_conditions(
        ax_slab, X, Y, T_slab, cmap_modified, x_extent, y_extent, phase_transition_depth_slab, bc_slab, "Slab Setup"
    )

    _ = draw_initial_conditions(ax_plume, X, Y, T_plume, cmap_modified, x_extent, y_extent, phase_transition_depth_plume, bc_plume, "Plume Setup")

    ticks = [-500, -250, 0, 250, 500]
    cax = fig.add_axes([0.30, -0.02, 0.55, 0.028])  # pyright: ignore
    cbar = fig.colorbar(contour_slab, cax=cax, ticks=ticks, orientation="horizontal")
    cbar.set_label("Thermal Anomaly (K)")
    cbar.ax.tick_params(labelsize=13)

    plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.05)


if __name__ == "__main__":
    main()
