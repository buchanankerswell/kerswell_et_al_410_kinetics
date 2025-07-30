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
    visualization_type = args.visualization_type

    if visualization_type == "2d_shell":
        active_plot_fields = [
            "nonadiabatic_temperature",
            "nonadiabatic_pressure",
            "nonadiabatic_density",
            "viscosity",
            "shear_stress_xx",
            "shear_stress_xy",
            "shear_stress_yy",
            "stress_second_invariant",
            "strain_rate",
            "seismic_Vp",
            "seismic_Vs",
        ]
    elif visualization_type == "2d_box":
        active_plot_fields = [
            # "T",
            # "adiabatic_temperature",
            "nonadiabatic_temperature",
            # "p",
            # "adiabatic_pressure",
            "nonadiabatic_pressure",
            "X_field",
            "density",
            # "adiabatic_density",
            "nonadiabatic_density",
            "velocity",
            # "stress_xx",
            # "stress_yy",
            # "shear_stress_xx",
            # "shear_stress_xy",
            # "shear_stress_yy",
            "stress_second_invariant",
            # "strain_rate",
            # "viscosity",
            "seismic_Vp",
            "seismic_Vs",
            # "Vp_anomaly",
            # "Vs_anomaly",
            "driving_force",
            "reaction_rate_C0",
        ]
    else:
        raise ValueError(f"Unrecognized visualization type {visualization_type}!")

    plot_config = PyVistaModelConfig()

    plot_config.file_mapping = {k: plot_config.file_mapping[k] for k in active_plot_fields if k in plot_config.file_mapping}

    plot_config.default_fig_dir = out_fig_dir

    if visualization_type == "2d_shell":
        plot_config.draw_mesh_plots = True
        plot_config.draw_centerline_depth_plots = False
        plot_config.draw_deflection_plots = False
        plot_config.camera_center_zoom = False

    elif visualization_type == "2d_box":
        plot_config.draw_mesh_plots = True
        plot_config.draw_centerline_depth_plots = True
        plot_config.draw_deflection_plots = True
        plot_config.camera_center_zoom = True
        plot_config.title_position = (0.39, 0.90)
        plot_config.cbar_position = [0.30, 0.05]
        plot_config.depth_contour_depths_km = []
        plot_config.depth_contour_tolerance_km = 1
        plot_config.depth_contour_line_widths = [3]
        plot_config.scale_bar_length_fraction = 0.253

        plot_config.clim_mapping = {
            "T": "auto",
            "adiabatic_temperature": "auto",
            "nonadiabatic_temperature": "auto",
            "p": "auto",
            "adiabatic_pressure": "auto",
            "nonadiabatic_pressure": "auto",
            "X_field": (0.0, 1.0),
            "density": "auto",
            "adiabatic_density": "auto",
            "nonadiabatic_density": "auto",
            "velocity": "auto",
            "stress_xx": "auto",
            "stress_yy": "auto",
            "shear_stress_xx": "auto",
            "shear_stress_xy": "auto",
            "shear_stress_yy": "auto",
            "stress_second_invariant": "auto",
            "strain_rate": "auto",
            "viscosity": "auto",
            "seismic_Vp": "auto",
            "seismic_Vs": "auto",
            "Vp_anomaly": "auto",
            "Vs_anomaly": "auto",
            "driving_force": "auto",
            "reaction_rate_C0": "auto",
        }
        plot_config.fmt_mapping = {
            "T": "%.0f",
            "adiabatic_temperature": "%.0f",
            "nonadiabatic_temperature": "%.0f",
            "p": "%.0f",
            "adiabatic_pressure": "%.0f",
            "nonadiabatic_pressure": "%.0f",
            "X_field": "%.1f",
            "density": "%.1f",
            "adiabatic_density": "%.2f",
            "nonadiabatic_density": "%.2f",
            "velocity": "%.1f",
            "stress_xx": "%.1f",
            "stress_yy": "%.1f",
            "shear_stress_xx": "%.0f",
            "shear_stress_xy": "%.0f",
            "shear_stress_yy": "%.0f",
            "principal_stress_2": "%.0f",
            "principal_stress_1": "%.0f",
            "stress_second_invariant": "%.0f",
            "strain_rate": "%.1f",
            "viscosity": "%.0f",
            "seismic_Vp": "%.1f",
            "seismic_Vs": "%.1f",
            "Vp_anomaly": "%.0f",
            "Vs_anomaly": "%.0f",
            "driving_force": "%.0f",
            "reaction_rate_C0": "%.1f",
        }
        plot_config.scale_mapping = {
            "T": 1,
            "adiabatic_temperature": 1,
            "nonadiabatic_temperature": 1,
            "p": 1e-9,
            "adiabatic_pressure": 1e-9,
            "nonadiabatic_pressure": 1e-6,
            "X_field": 1,
            "density": 1e-3,
            "adiabatic_density": 1e-3,
            "nonadiabatic_density": 1e-3,
            "velocity": 1e2,
            "stress_xx": 1e-9,
            "stress_yy": 1e-9,
            "shear_stress_xx": 1e-6,
            "shear_stress_xy": 1e-6,
            "shear_stress_yy": 1e-6,
            "principal_stress_1": 1e-6,
            "principal_stress_2": 1e-6,
            "stress_second_invariant": 1e-6,
            "strain_rate": 1,
            "viscosity": 1,
            "seismic_Vp": 1e-3,
            "seismic_Vs": 1e-3,
            "Vp_anomaly": 1,
            "Vs_anomaly": 1,
            "driving_force": 1e-3,
            "reaction_rate_C0": 3.154e13,
        }

        plot_config.bar_mapping = {
            "T": "$T$ (K)",
            "adiabatic_temperature": "$\\bar{T}$ (K)",
            "nonadiabatic_temperature": "$\\hat{T}$ (K)",
            "p": "$P$ (GPa)",
            "adiabatic_pressure": "$\\bar{P}$ (GPa)",
            "nonadiabatic_pressure": "$\\hat{P}$ (MPa)",
            "X_field": "X",
            "density": "$\\rho$ (g/cm$^3$)",
            "adiabatic_density": "$\\bar{\\rho}$ (g/cm$^3$)",
            "nonadiabatic_density": "$\\hat{\\rho}$ (g/cm$^3$)",
            "velocity": "$\\vec{u}$ (cm/yr)",
            "stress_xx": "$\\sigma_{xx}$ (GPa)",
            "stress_yy": "$\\sigma_{yy}$ (GPa)",
            "shear_stress_xx": "$\\sigma^{\\prime}_{xx}$ (MPa)",
            "shear_stress_xy": "$\\sigma^{\\prime}_{xy}$ (MPa)",
            "shear_stress_yy": "$\\sigma^{\\prime}_{yy}$ (MPa)",
            "principal_stress_1": "$\\sigma^{\\prime}_1$ (MPa)",
            "principal_stress_2": "$\\sigma^{\\prime}_2$ (MPa)",
            "stress_second_invariant": "$\\sigma_{II}$ (MPa)",
            "strain_rate": "Log $\\dot{\\epsilon}_{II}$ (1/s)",
            "viscosity": "Log $\\eta$ (Pa s)",
            "seismic_Vp": "$V_p$ (km/s)",
            "seismic_Vs": "$V_s$ (km/s)",
            "Vp_anomaly": "$\\hat{V}_p$ (%)",
            "Vs_anomaly": "$\\hat{V}_s$ (%)",
            "driving_force": "$\\Delta G$ (kJ/mol)",
            "reaction_rate_C0": "$dX/dt$ (1/Ma)",
        }

        plot_config.velocity_mapping = {
            "T": False,
            "adiabatic_temperature": False,
            "nonadiabatic_temperature": True,
            "p": False,
            "adiabatic_pressure": False,
            "nonadiabatic_pressure": False,
            "X_field": False,
            "density": False,
            "adiabatic_density": False,
            "nonadiabatic_density": False,
            "velocity": True,
            "stress_xx": False,
            "stress_yy": False,
            "shear_stress_xx": False,
            "shear_stress_xy": False,
            "shear_stress_yy": False,
            "principal_stress_1": False,
            "principal_stress_2": False,
            "stress_second_invariant": False,
            "strain_rate": False,
            "viscosity": False,
            "seismic_Vp": False,
            "seismic_Vs": False,
            "Vp_anomaly": False,
            "Vs_anomaly": False,
            "driving_force": False,
            "reaction_rate_C0": False,
        }
    else:
        raise ValueError(f"Unrecognized visualiztion type {visualization_type}!")

    pvtu_in_dirs = dict(zip(model_ids, in_dirs))

    visualizer = PyVistaModelVisualizer(
        plot_config=plot_config,
        pvtu_in_dirs=pvtu_in_dirs,
        tsteps_mesh=tsteps_mesh,
        tsteps_profile=tsteps_profile,
        tsteps_topography=tsteps_topography,
    )

    visualizer.draw()

    tile_sets = {
        "set0": {
            "tags": None,
            "fields": [
                "nonadiabatic_temperature",
                "reaction_rate_C0",
                "X_field",
            ],
        },
        "set1": {
            "tags": None,
            "fields": [
                "nonadiabatic_temperature",
                "driving_force",
                "reaction_rate_C0",
            ],
        },
        "set2": {
            "tags": None,
            "fields": [
                "nonadiabatic_temperature",
                "X_field",
                "nonadiabatic_density",
            ],
        },
        "set3": {
            "tags": None,
            "fields": [
                "velocity",
                "nonadiabatic_pressure",
                "stress_second_invariant",
            ],
        },
        "set4": {
            "tags": None,
            "fields": [
                "density",
                "seismic_Vp",
                "seismic_Vs",
            ],
        },
    }

    print("    --------------------------------------------------")
    print("    Tiling images")
    print("    --------------------------------------------------")

    for model_id in pvtu_in_dirs.keys():
        for config in tile_sets.values():
            fields = config.get("fields", None)
            tags = config.get("tags", None)

            if fields and len(fields) == 3:
                tiler = ImageTiler(
                    plot_config=plot_config,
                    out_fig_dir=out_fig_dir / model_id,
                    field1=fields[0],
                    field2=fields[1],
                    field3=fields[2],
                    movie=False,
                    fps=12,
                    tags=tags,
                )
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
    parser.add_argument("--visualization-type", type=str, help="Visualization type")

    return parser.parse_args()


if __name__ == "__main__":
    main()
