#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
from argparse import ArgumentParser, Namespace
from pathlib import Path

from tile import ImageTiler
from visualize import PyVistaModelConfig, PyVistaModelVisualizer


#######################################################
## .1. Main                                      !!! ##
#######################################################
def main():
    """"""
    args = parse_arguments()

    model_ids = args.model_ids or []
    tsteps_mesh = args.timesteps_mesh or []
    tsteps_profile = args.timesteps_profile or []
    tsteps_topography = args.timesteps_topography or []
    in_dirs = [Path(p) for p in args.in_dirs] if args.in_dirs else []
    out_fig_dir = Path(args.out_fig_dir) if args.out_fig_dir else Path("./figures")

    active_plot_fields = [
        "nonadiabatic_temperature",
        "nonadiabatic_pressure",
        "nonadiabatic_density",
        "viscosity",
        "stress_second_invariant",
        "strain_rate",
        "velocity",
        "seismic_Vp",
        "seismic_Vs",
        "arrhenius",
        "thermodynamic",
        "reaction_rate_C0",
        "X_field",
    ]

    plot_config = PyVistaModelConfig()
    plot_config.file_mapping = {k: plot_config.file_mapping[k] for k in active_plot_fields if k in plot_config.file_mapping}
    plot_config.default_fig_dir = out_fig_dir
    pvtu_in_dirs = dict(zip(model_ids, in_dirs))

    visualizer = PyVistaModelVisualizer(plot_config, pvtu_in_dirs, tsteps_mesh, tsteps_profile, tsteps_topography)
    visualizer.draw()

    tile_sets = {
        "set0": {"tags": "abc", "fields": ["nonadiabatic_temperature", "nonadiabatic_pressure", "nonadiabatic_density"]},
        "set1": {"tags": "def", "fields": ["arrhenius", "thermodynamic", "reaction_rate_C0", ]},
        "set3": {"tags": "ghi", "fields": ["X_field", "seismic_Vp", "seismic_Vs", ]},
        "set4": {"tags": "abc", "fields": ["nonadiabatic_temperature", "reaction_rate_C0", "X_field"]},
        "set5": {"tags": "def", "fields": ["nonadiabatic_temperature", "reaction_rate_C0", "X_field"]},
        "set6": {"tags": "ghi", "fields": ["nonadiabatic_temperature", "reaction_rate_C0", "X_field"]},
        "set7": {"tags": "abc", "fields": ["nonadiabatic_temperature", "nonadiabatic_density", "seismic_Vp"]},
        "set8": {"tags": "def", "fields": ["nonadiabatic_temperature", "nonadiabatic_density", "seismic_Vp"]},
        "set9": {"tags": "ghi", "fields": ["nonadiabatic_temperature", "nonadiabatic_density", "seismic_Vp"]},
        "set10": {"tags": "jkl", "fields": ["viscosity", "stress_second_invariant", "strain_rate"]},
    }

    print("    --------------------------------------------------")
    print("    Tiling images")
    print("    --------------------------------------------------")

    for model_id in pvtu_in_dirs.keys():
        for config in tile_sets.values():
            fields = config.get("fields", None)
            tags = config.get("tags", None)

            if fields and len(fields) == 3:
                tiler = ImageTiler(plot_config, out_fig_dir / "meshes" / model_id, fields[0], fields[1], fields[2], tags)
                tiler.tile_images()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_arguments() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Visualize simulation results.")
    parser.add_argument("--model-ids", nargs="+", type=str, help="Simulation names")
    parser.add_argument("--timesteps-mesh", nargs="+", type=int, help="Timesteps to visualize mesh")
    parser.add_argument("--timesteps-profile", nargs="+", type=int, help="Timesteps to visualize profile")
    parser.add_argument("--timesteps-topography", nargs="+", type=int, help="Timesteps to visualize topography")
    parser.add_argument("--in-dirs", nargs="+", type=str, help="Simulation results directories")
    parser.add_argument("--out-fig-dir", type=str, help="Output figure directory")

    return parser.parse_args()


if __name__ == "__main__":
    main()
