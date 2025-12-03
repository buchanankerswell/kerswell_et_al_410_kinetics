#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
import gc
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv
from PIL import Image
from scipy.signal import savgol_filter
from scipy.spatial import ConvexHull

__import__("vtk")


#######################################################
## .1. PyVistaModelConfig                        !!! ##
#######################################################
@dataclass
class PyVistaModelConfig:
    """Holds all configuration for PyVista visualizations."""

    draw_mesh_plots: bool = True
    draw_depth_profile_plots: bool = False
    draw_topography_plots: bool = False
    draw_X_field_contours: bool = True

    # Mappings
    file_mapping: dict[str, str] = field(default_factory=dict)
    cmap_mapping: dict[str, str] = field(default_factory=dict)
    clim_mapping: dict[str, tuple[float, float] | str] = field(default_factory=dict)
    bar_mapping: dict[str, str] = field(default_factory=dict)
    fmt_mapping: dict[str, str] = field(default_factory=dict)
    scale_mapping: dict[str, float] = field(default_factory=dict)
    velocity_mapping: dict[str, bool] = field(default_factory=dict)

    stress_components: dict[str, int] = field(
        default_factory=lambda: {
            "stress_xx": 0,
            "stress_xy": 1,
            "stress_xz": 2,
            "stress_yx": 3,
            "stress_yy": 4,
            "stress_yz": 5,
            "stress_zx": 6,
            "stress_zy": 7,
            "stress_zz": 8,
        }
    )
    shear_stress_components: dict[str, int] = field(
        default_factory=lambda: {
            "shear_stress_xx": 0,
            "shear_stress_xy": 1,
            "shear_stress_xz": 2,
            "shear_stress_yx": 3,
            "shear_stress_yy": 4,
            "shear_stress_yz": 5,
            "shear_stress_zx": 6,
            "shear_stress_zy": 7,
            "shear_stress_zz": 8,
        }
    )
    add_mesh_fields: list[str] = field(
        default_factory=lambda: [
            "differential_stress",
            "nonadiabatic_density",
        ]
    )

    n_colors: int = 11
    show_edges: bool = False
    default_fig_dir: Path = Path("figs/simulation")
    plotter_window_size: list[int] = field(default_factory=lambda: [1920, 1920])
    plotter_background: str = "#FFFFFF"
    plotter_lighting: str = "none"
    edge_opacity: float = 0.3
    cbar_vertical: bool = False
    cbar_title_font_size: int = 115
    cbar_label_font_size: int = 115
    cbar_width: float = 0.6
    cbar_height: float = 0.13
    cbar_n_labels: int = 3
    cbar_position: list[float] = field(default_factory=lambda: [0.20, 0.12])
    camera_center_zoom: bool = True
    camera_full_view: bool = True
    camera_y_shift_factor: float = -0.15
    camera_zoom_factor: float = 1.95
    screenshot_dpi: tuple[int, int] = field(default_factory=lambda: (330, 330))
    velocity_glyph_factor: float = 18e5
    stress_glyph_factor: float = 1e-3
    glyph_line_width: float = 4
    scale_bar_enabled: bool = True
    scale_bar_color: str = "black"
    scale_bar_thickness: float = 30
    scale_bar_length_fraction: float = 0.222
    scale_bar_position: tuple = (0.05, 0.08)
    scale_bar_label_font_size: int = 115
    scale_bar_label_shift_factor: float = 0.01

    topography_plot_width: float = 7.0
    topography_plot_height: float = 5.5

    depth_profile_plot_width: float = 7.0
    depth_profile_plot_height: float = 5.5
    depth_profile_surface_depth: float = 280e3
    depth_profile_reaction_depth: float = 410e3
    depth_profile_plot_rcParams: dict[str, Any] = field(
        default_factory=lambda: {
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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __post_init__(self):
        if not self.file_mapping:
            self.file_mapping = {
                "T": "temperature-full",
                "adiabatic_temperature": "temperature-adiabatic",
                "nonadiabatic_temperature": "temperature-nonadiabatic",
                "p": "pressure-full",
                "nonadiabatic_pressure": "pressure-nonadiabatic",
                "X_field": "X-field",
                "density": "density-full",
                "nonadiabatic_density": "density-nonadiabatic",
                "velocity": "velocity-cart",
                "stress_xx": "sigma-full-xx",
                "stress_yy": "sigma-full-yy",
                "shear_stress_xx": "sigma-deviatoric-xx",
                "shear_stress_xy": "sigma-deviatoric-xy",
                "shear_stress_yy": "sigma-deviatoric-yy",
                "stress_second_invariant": "sigma-ii",
                "strain_rate": "strain-rate",
                "adiabatic_pressure": "pressure-adiabatic",
                "adiabatic_density": "density-adiabatic",
                "adiabatic_density_derivative": "density-adiabatic-gradient",
                "specific_heat": "specific-heat-capacity",
                "thermal_expansivity": "thermal-expansivity",
                "compressibility": "compressibility",
                "viscosity": "viscosity",
                "principal_stress_1": "sigma-one",
                "principal_stress_2": "sigma-two",
                "differential_stress": "sigma-delta",
                "seismic_Vp": "vp",
                "seismic_Vs": "vs",
                "arrhenius": "arrhenius",
                "thermodynamic": "thermodynamic",
                "reaction_rate_C0": "reaction-rate",
            }
        if not self.cmap_mapping:
            self.cmap_mapping = {
                "T": "magma",
                "adiabatic_temperature": "magma",
                "nonadiabatic_temperature": "seismic",
                "p": "Oranges",
                "adiabatic_pressure": "Oranges",
                "nonadiabatic_pressure": "PuOr_r",
                "X_field": "GnBu",
                "density": "Purples",
                "adiabatic_density": "Purples",
                "nonadiabatic_density": "BrBG",
                "adiabatic_density_derivative": "Purples",
                "specific_heat": "bone_r",
                "thermal_expansivity": "pink_r",
                "compressibility": "gist_heat_r",
                "viscosity": "viridis_r",
                "velocity": "BrBG",
                "stress_xx": "BuPu",
                "stress_yy": "BuPu",
                "shear_stress_xx": "PiYG",
                "shear_stress_xy": "PiYG",
                "shear_stress_yy": "PiYG",
                "principal_stress_1": "PRGn",
                "principal_stress_2": "PRGn",
                "differential_stress": "Oranges",
                "stress_second_invariant": "Reds",
                "strain_rate": "BuPu_r",
                "seismic_Vp": "copper_r",
                "seismic_Vs": "copper_r",
                "arrhenius": "bone_r",
                "thermodynamic": "pink_r",
                "reaction_rate_C0": "RdGy_r",
            }
        if not self.clim_mapping:
            self.clim_mapping = {
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
                "adiabatic_density_derivative": "auto",
                "specific_heat": "auto",
                "thermal_expansivity": "auto",
                "compressibility": "auto",
                "viscosity": "auto",
                "velocity": "auto",
                "stress_xx": "auto",
                "stress_yy": "auto",
                "shear_stress_xx": "auto",
                "shear_stress_xy": "auto",
                "shear_stress_yy": "auto",
                "principal_stress_1": "auto",
                "principal_stress_2": "auto",
                "differential_stress": "auto",
                "stress_second_invariant": "auto",
                "strain_rate": "auto",
                "seismic_Vp": "auto",
                "seismic_Vs": "auto",
                "arrhenius": "auto",
                "thermodynamic": "auto",
                "reaction_rate_C0": "auto",
            }
        if not self.bar_mapping:
            self.bar_mapping = {
                "T": "$T$ (K)",
                "adiabatic_temperature": "$\\bar{T}$ (K)",
                "nonadiabatic_temperature": "$\\hat{T}$ (K)",
                "p": "$P$ (GPa)",
                "adiabatic_pressure": "$\\bar{P}$ (MPa)",
                "nonadiabatic_pressure": "$\\hat{P}$ (MPa)",
                "X_field": "X",
                "density": "$\\rho$ (g/cm$^3$)",
                "adiabatic_density": "$\\bar{\\rho}$ (g/cm$^3$)",
                "nonadiabatic_density": "$\\hat{\\rho}$ (g/cm$^3$)",
                "adiabatic_density_derivative": "$\\nabla \\rho$ (kg/m$^3$ km)",
                "specific_heat": "$C_p$ (kJ/K kg)",
                "thermal_expansivity": "$\\alpha$ (1/K)",
                "compressibility": "$\\beta_S$ (1/Pa)",
                "viscosity": "Log $\\eta$ (Pa s)",
                "velocity": "$\\vec{u}$ (cm/yr)",
                "stress_xx": "$\\sigma_{xx}$ (GPa)",
                "stress_yy": "$\\sigma_{yy}$ (GPa)",
                "shear_stress_xx": "$\\sigma^{\\prime}_{xx}$ (GPa)",
                "shear_stress_xy": "$\\sigma^{\\prime}_{xy}$ (GPa)",
                "shear_stress_yy": "$\\sigma^{\\prime}_{yy}$ (GPa)",
                "principal_stress_1": "$\\sigma^{\\prime}_1$ (GPa)",
                "principal_stress_2": "$\\sigma^{\\prime}_2$ (GPa)",
                "differential_stress": "$\\sigma^{\\prime}_1$-$\\sigma^{\\prime}_2$ (GPa)",
                "stress_second_invariant": "$\\sigma_{II}$ (MPa)",
                "strain_rate": "Log $\\dot{\\epsilon}_{II}$ (1/s)",
                "seismic_Vp": "$V_p$ (km/s)",
                "seismic_Vs": "$V_s$ (km/s)",
                "arrhenius": "Log Arrhenius term",
                "thermodynamic": "Thermodynamic term",
                "reaction_rate_C0": "Log $\\dot{X}$ (1/Ma)",
            }
        if not self.fmt_mapping:
            self.fmt_mapping = {
                "T": "%.0f",
                "adiabatic_temperature": "%.0f",
                "nonadiabatic_temperature": "%.0f",
                "p": "%.0f",
                "adiabatic_pressure": "%.0f",
                "nonadiabatic_pressure": "%.0f",
                "X_field": "%.1f",
                "density": "%.1f",
                "adiabatic_density": "%.1f",
                "nonadiabatic_density": "%.2f",
                "adiabatic_density_derivative": "%.0f",
                "specific_heat": "%.2f",
                "thermal_expansivity": "%.1f",
                "compressibility": "%.1e",
                "viscosity": "%.1f",
                "velocity": "%.1f",
                "stress_xx": "%.0f",
                "stress_yy": "%.0f",
                "shear_stress_xx": "%.1f",
                "shear_stress_xy": "%.1f",
                "shear_stress_yx": "%.1f",
                "shear_stress_yy": "%.1f",
                "principal_stress_2": "%.1f",
                "principal_stress_1": "%.1f",
                "differential_stress": "%.1f",
                "stress_second_invariant": "%.0f",
                "strain_rate": "%.1f",
                "seismic_Vp": "%.1f",
                "seismic_Vs": "%.1f",
                "arrhenius": "%.1f",
                "thermodynamic": "%.1f",
                "reaction_rate_C0": "%.2f",
            }
        if not self.scale_mapping:
            self.scale_mapping = {
                "T": 1,
                "adiabatic_temperature": 1,
                "nonadiabatic_temperature": 1,
                "p": 1e-9,
                "adiabatic_pressure": 1e-6,
                "nonadiabatic_pressure": 1e-6,
                "X_field": 1,
                "density": 1e-3,
                "adiabatic_density": 1e-3,
                "nonadiabatic_density": 1e-3,
                "adiabatic_density_derivative": 1e3,
                "specific_heat": 1e-3,
                "thermal_expansivity": 1e5,
                "compressibility": 1,
                "viscosity": 1,
                "velocity": 1e2,
                "stress_xx": 1e-9,
                "stress_yy": 1e-9,
                "shear_stress_xx": 1e-9,
                "shear_stress_xy": 1e-9,
                "shear_stress_yy": 1e-9,
                "principal_stress_1": 1e-9,
                "principal_stress_2": 1e-9,
                "differential_stress": 1e-9,
                "stress_second_invariant": 1e-6,
                "strain_rate": 3.154e13,
                "seismic_Vp": 1e-3,
                "seismic_Vs": 1e-3,
                "arrhenius": 1,
                "thermodynamic": 1,
                "reaction_rate_C0": 3.154e13,
            }
        if not self.velocity_mapping:
            self.velocity_mapping = {
                "T": False,
                "adiabatic_temperature": False,
                "nonadiabatic_temperature": True,
                "p": False,
                "nonadiabatic_pressure": False,
                "X_field": False,
                "density": False,
                "nonadiabatic_density": False,
                "velocity": False,
                "stress_xx": False,
                "stress_yy": False,
                "shear_stress_xx": False,
                "shear_stress_xy": False,
                "shear_stress_yy": False,
                "stress_second_invariant": False,
                "strain_rate": False,
                "adiabatic_pressure": False,
                "adiabatic_density": False,
                "adiabatic_density_derivative": False,
                "specific_heat": False,
                "thermal_expansivity": False,
                "compressibility": False,
                "viscosity": False,
                "principal_stress_1": False,
                "principal_stress_2": False,
                "differential_stress": False,
                "seismic_Vp": False,
                "seismic_Vs": False,
                "arrhenius": False,
                "thermodynamic": False,
                "reaction_rate_C0": False,
            }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def get_extra_fields(self) -> list[str]:
        """Returns a list of fields that can be derived from other mesh fields."""
        return list(self.stress_components.keys()) + list(self.shear_stress_components.keys()) + self.add_mesh_fields


#######################################################
## .2. PyVistaModelVisualizer                    !!! ##
#######################################################
@dataclass
class PyVistaModelVisualizer:
    plot_config: PyVistaModelConfig
    pvtu_in_dirs: dict[str, Path] = field(default_factory=dict)
    tsteps_mesh: list[int] = field(default_factory=list)
    tsteps_profile: list[int] = field(default_factory=list)
    tsteps_topography: list[int] = field(default_factory=list)
    verbosity: int = 0

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def draw(self):
        """Main loop to process all models and generate visualizations."""
        warnings.filterwarnings("ignore", category=RuntimeWarning, module="pyvista")

        cfg = self.plot_config
        model_out_data = self._get_ordered_pvtu_files()

        depth_profile_cache = {}
        depth_profile_summary = []
        reaction_depth = cfg.depth_profile_reaction_depth
        lower_bound = np.nan
        upper_bound = np.nan

        for model_id, (pvtu_files, timesteps) in model_out_data.items():
            print("    --------------------------------------------------")
            print(f"    Drawing mesh for: {model_id}")
            print("    --------------------------------------------------")

            if not pvtu_files:
                if self.verbosity >= 1:
                    print(f" !! Warning: no result files to process for model: {model_id}")
                continue

            out_fig_dir_mesh = cfg.default_fig_dir / "meshes" / model_id
            if not out_fig_dir_mesh.exists():
                out_fig_dir_mesh.mkdir(parents=True, exist_ok=True)

            Z_factor, B_factor = self._extract_z_b_values(model_id)

            for pvtu_path, tstep in zip(pvtu_files, timesteps):
                if tstep == -1 or tstep not in set(self.tsteps_mesh + self.tsteps_profile + self.tsteps_topography):
                    continue

                try:
                    mesh = pv.UnstructuredGrid(pvtu_path)
                    if mesh is None:
                        if self.verbosity >= 1:
                            print(f" !! Warning: pyvista.read returned None for pvtu_file:\n    {pvtu_path.name}")
                        continue
                except Exception as e:
                    if self.verbosity >= 1:
                        print(f" !! Warning: failed to read mesh pvtu_file {pvtu_path.name}:\n    {e}")
                    continue

                for field_name in cfg.file_mapping.keys():
                    self._prepare_mesh(mesh, field_name)
                    if not self._check_mesh(mesh, field_name):
                        return

                    field_id = cfg.file_mapping.get(field_name, field_name.replace("_", "-"))
                    time_myr = self._get_mesh_time_myr(mesh, model_id, tstep)
                    cmap, clim = self._configure_cmap(mesh, field_name, 0.0)
                    x_pos = np.nan

                    if tstep in self.tsteps_profile or tstep in self.tsteps_topography:
                        depth_profile_cache.setdefault(model_id, {}).setdefault(tstep, {})
                        depth_profile_cache[model_id][tstep]["time_myr"] = time_myr
                        depth_profile_cache[model_id][tstep].setdefault("fields", {})

                        if "X_field" not in mesh.point_data:
                            continue

                        displacement, width, x_pos = self._find_vertical_profile_for_410_structure(
                            mesh.copy(), "X_field", reaction_depth, Z_factor, B_factor
                        )
                        depth_profile_cache[model_id][tstep]["displacement"] = displacement
                        depth_profile_cache[model_id][tstep]["width"] = width
                        depth_profile_cache[model_id][tstep]["profile_x"] = x_pos

                        if field_name in mesh.point_data:
                            depths, values = self._get_depth_profile_data_at_x_pos(mesh.copy(), field_name, x_pos)
                            depth_profile_cache[model_id][tstep]["fields"][field_name] = {"depths": depths, "values": values, "clim": clim}

                            if field_name in ["reaction_rate_C0", "velocity"]:
                                key = "max_reaction_rate" if field_name == "reaction_rate_C0" else "max_velocity"

                                if width <= 5e3:
                                    lower_factor = 2
                                    upper_factor = 1
                                else:
                                    lower_factor = 1.1
                                    upper_factor = 0.1

                                lower_bound = max(np.nanmin(depths), reaction_depth + displacement - width * lower_factor)
                                upper_bound = min(np.nanmax(depths), reaction_depth + displacement + width * upper_factor)
                                mask_depth = (depths >= lower_bound) & (depths <= upper_bound)

                                if np.any(mask_depth):
                                    depth_profile_cache[model_id][tstep][key] = np.nanmax(abs(values[mask_depth]))
                                else:
                                    depth_profile_cache[model_id][tstep][key] = np.nan

                        if all(k in depth_profile_cache[model_id][tstep] for k in ["displacement", "max_reaction_rate", "max_velocity"]):
                            depth_profile_summary.append(
                                {
                                    "model_id": model_id,
                                    "Z_factor": Z_factor,
                                    "B_factor": B_factor,
                                    "timestep": tstep,
                                    "time_myr": time_myr,
                                    "displacement": depth_profile_cache[model_id][tstep]["displacement"],
                                    "width": depth_profile_cache[model_id][tstep]["width"],
                                    "max_velocity": depth_profile_cache[model_id][tstep]["max_velocity"],
                                    "max_reaction_rate": depth_profile_cache[model_id][tstep]["max_reaction_rate"],
                                    "profile_x": depth_profile_cache[model_id][tstep]["profile_x"],
                                }
                            )

                    if cfg.draw_mesh_plots and tstep in self.tsteps_mesh:
                        if field_name == "X_field":
                            draw_vertical_profile = True
                            draw_horizontal_window = True
                            lb = lower_bound - cfg.depth_profile_surface_depth if not np.isnan(lower_bound) else np.nan
                            ub = upper_bound - cfg.depth_profile_surface_depth if not np.isnan(upper_bound) else np.nan
                        else:
                            draw_vertical_profile = False
                            draw_horizontal_window = False
                            lb = np.nan
                            ub = np.nan

                        out_path = out_fig_dir_mesh / f"{model_id.replace('_', '-')}-{field_id}-{str(tstep).zfill(4)}.png"
                        self.draw_mesh(mesh.copy(), field_name, cmap, clim, x_pos, draw_vertical_profile, draw_horizontal_window, lb, ub, out_path)

        if depth_profile_summary:
            try:
                df = pd.DataFrame(depth_profile_summary)
                df = df.drop_duplicates()
                csv_path = Path("simulation") / "data" / "depth-profile-data.csv"

                if not csv_path.exists():
                    csv_path.parent.mkdir(parents=True, exist_ok=True)
                    df.to_csv(csv_path, index=False)
            except Exception as e:
                print(f"!! Error: Failed to write CSV file with pandas:\n    {e}")

        if cfg.draw_depth_profile_plots and depth_profile_cache:
            print("    --------------------------------------------------")
            print("    Drawing 410 displacements")
            print("    --------------------------------------------------")

            for model_type in ["plume", "slab"]:
                selected_models = [m for m in depth_profile_cache.keys() if model_type in m]
                selected_models = sorted(selected_models, key=self._extract_sort_key)

                displacement_data_grouped = {}
                for model_id, tsteps in depth_profile_cache.items():
                    if model_id in selected_models:
                        times = []
                        defs = []
                        sorted_tsteps = sorted(tsteps.keys())
                        for tstep in sorted_tsteps:
                            if tstep in self.tsteps_topography:
                                data = tsteps[tstep]
                                if "displacement" in data and "time_myr" in data:
                                    times.append(data["time_myr"])
                                    defs.append(data["displacement"])
                        if times and defs:
                            displacement_data_grouped[model_id] = {
                                "time_myr": np.array(times, dtype=np.float32),
                                "displacement": np.array(defs, dtype=np.float32),
                            }

                out_fig_dir_displacement = cfg.default_fig_dir / "topography"
                out_fig_dir_displacement.mkdir(parents=True, exist_ok=True)
                out_path = out_fig_dir_displacement / f"{model_type}-displacement.png"
                self.draw_displacement(displacement_data_grouped, model_type, out_path)

        if cfg.draw_depth_profile_plots and depth_profile_cache:
            print("    --------------------------------------------------")
            print("    Drawing 410 widths")
            print("    --------------------------------------------------")

            for model_type in ["plume", "slab"]:
                selected_models = [m for m in depth_profile_cache.keys() if model_type in m]
                selected_models = sorted(selected_models, key=self._extract_sort_key)

                width_data_grouped = {}
                for model_id, tsteps in depth_profile_cache.items():
                    if model_id in selected_models:
                        times = []
                        defs = []
                        sorted_tsteps = sorted(tsteps.keys())
                        for tstep in sorted_tsteps:
                            if tstep in self.tsteps_topography:
                                data = tsteps[tstep]
                                if "width" in data and "time_myr" in data:
                                    times.append(data["time_myr"])
                                    defs.append(data["width"])
                        if times and defs:
                            width_data_grouped[model_id] = {
                                "time_myr": np.array(times, dtype=np.float32),
                                "width": np.array(defs, dtype=np.float32),
                            }

                out_fig_dir_width = cfg.default_fig_dir / "topography"
                out_fig_dir_width.mkdir(parents=True, exist_ok=True)
                out_path = out_fig_dir_width / f"{model_type}-width.png"
                self.draw_width(width_data_grouped, model_type, out_path)

        if cfg.draw_depth_profile_plots and depth_profile_cache:
            print("    --------------------------------------------------")
            print("    Drawing depth profiles")
            print("    --------------------------------------------------")

            for model_type in ["plume", "slab"]:
                selected_models = [m for m in depth_profile_cache.keys() if model_type in m]
                selected_models = sorted(selected_models, key=self._extract_sort_key)

                all_fields = set()
                for m in selected_models:
                    for tstep in depth_profile_cache[m].keys():
                        if tstep in self.tsteps_profile:
                            all_fields.update(depth_profile_cache[m][tstep].get("fields", {}).keys())

                for field_name in sorted(all_fields):
                    all_timesteps = sorted(
                        {t for m in selected_models for t in depth_profile_cache[m] if field_name in depth_profile_cache[m][t].get("fields", {})}
                    )

                    for tstep in all_timesteps:
                        if tstep in self.tsteps_profile:
                            profile_data_grouped = {}

                            time_myr = 0
                            for model_id in selected_models:
                                if tstep not in depth_profile_cache[model_id]:
                                    continue
                                field_data = depth_profile_cache[model_id][tstep].get("fields", {}).get(field_name)
                                if not field_data:
                                    continue
                                profile_data_grouped[model_id] = field_data
                                time_myr = depth_profile_cache[model_id][tstep]["time_myr"]

                            if not profile_data_grouped:
                                continue

                            out_fig_dir_profiles = cfg.default_fig_dir / "profiles"
                            out_fig_dir_profiles.mkdir(parents=True, exist_ok=True)

                            field_id = cfg.file_mapping.get(field_name, field_name.replace("_", "-"))
                            type_id = f"{model_type}"
                            out_path = out_fig_dir_profiles / f"{type_id}-{field_id}-depth-profile-{str(tstep).zfill(4)}.png"

                            self.draw_depth_profile(profile_data_grouped, field_name, model_type, out_path, reaction_depth)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def draw_mesh(
        self,
        mesh: pv.UnstructuredGrid,
        field_name: str,
        cmap: mcolors.Colormap | mcolors.ListedColormap,
        clim: tuple[float, float],
        x_pos: float,
        draw_vertical_profile: bool,
        draw_horizontal_window: bool,
        lower_bound: float,
        upper_bound: float,
        out_path: Path,
    ) -> None:
        """Visualizes a single field on the mesh using PyVista."""
        if out_path.exists():
            print(f" -- Found plot: {out_path.name}!")
            return

        print(f" -> {out_path.name}")

        cfg = self.plot_config

        plotter: pv.Plotter = pv.Plotter(
            off_screen=True,
            window_size=cfg.plotter_window_size,
            lighting=cfg.plotter_lighting,  # type: ignore[assignment]
        )
        plotter.ren_win.OffScreenRenderingOn()
        plotter.set_background(cfg.plotter_background)  # type: ignore

        glyph_arrows = None
        if (field_name == "velocity" or cfg.velocity_mapping[field_name]) and "velocity" in mesh.point_data:
            if "velocity_mag" not in mesh.point_data:
                mesh["velocity_mag"] = np.linalg.norm(mesh["velocity"], axis=1)
            glyph_arrows = mesh.glyph(
                orient="velocity",  # type: ignore[arg-type]
                scale="velocity_mag",  # type: ignore[arg-type]
                factor=cfg.velocity_glyph_factor,
                geom=pv.Arrow(),
                tolerance=0.05,
            )
        elif field_name == "principal_stress_1" and "principal_stress_direction_1" in mesh.point_data:
            glyph_arrows = mesh.glyph(
                orient="principal_stress_direction_1",  # type: ignore[arg-type]
                scale="principal_stress_1",  # type: ignore[arg-type]
                factor=cfg.stress_glyph_factor,
                geom=pv.Line(),
            )
        elif field_name == "principal_stress_2" and "principal_stress_direction_2" in mesh.point_data:
            glyph_arrows = mesh.glyph(
                orient="principal_stress_direction_2",  # type: ignore[arg-type]
                scale="principal_stress_2",  # type: ignore[arg-type]
                factor=cfg.stress_glyph_factor,
                geom=pv.Line(),
            )

        if glyph_arrows:
            plotter.add_mesh(
                glyph_arrows,
                color="black",
                line_width=cfg.glyph_line_width,
                render_lines_as_tubes=False,
            )

        sargs = dict(
            title=cfg.bar_mapping.get(field_name, field_name),
            vertical=cfg.cbar_vertical,
            title_font_size=cfg.cbar_title_font_size,
            label_font_size=cfg.cbar_label_font_size,
            fmt=cfg.fmt_mapping.get(field_name, "%.1f"),
            width=cfg.cbar_width,
            height=cfg.cbar_height,
            n_labels=cfg.cbar_n_labels,
            position_x=cfg.cbar_position[0],
            position_y=cfg.cbar_position[1],
        )

        plotter.add_mesh(
            mesh,
            scalars=f"{field_name}_viz",
            cmap=cmap,
            clim=clim,
            show_edges=cfg.show_edges,
            edge_opacity=cfg.edge_opacity,
            scalar_bar_args=sargs,  # type: ignore[assignment]
            nan_color="#FEFEFE",
        )

        if self.plot_config.draw_X_field_contours:
            self._add_X_field_contours(plotter, mesh)

        camera_settings = self._compute_camera_settings(mesh, cfg.camera_full_view)
        plotter.camera_position = camera_settings
        plotter.enable_parallel_projection()

        if cfg.scale_bar_enabled:
            self._add_scale_bar(mesh, plotter)

        if not np.isnan(lower_bound) and not np.isnan(upper_bound) and draw_horizontal_window:
            self._add_depth_contour(plotter, mesh, lower_bound)
            self._add_depth_contour(plotter, mesh, upper_bound)

        if not np.isnan(x_pos) and draw_vertical_profile:
            self._add_vertical_contour(plotter, mesh, x_pos)

        gc.collect()
        plotter.screenshot(out_path)
        gc.collect()
        plotter.close()
        pv.close_all()

        del mesh, plotter, glyph_arrows
        gc.collect()

        try:
            img = Image.open(out_path)
            img.save(out_path, dpi=cfg.screenshot_dpi)
        except Exception as e:
            print(f" !! Error: PIL failed to resave {out_path} with new DPI:\n    {e}")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def draw_depth_profile(
        self,
        profile_data_grouped: dict[str, dict],
        field_name: str,
        model_type: str,
        out_path: Path,
        reaction_depth: float,
        y_offset: float = 150e3,
        depth_points: int = 200,
        smoothing_range: float = 2e3,
    ) -> None:
        """"""
        if out_path.exists():
            print(f" -- Found plot: {out_path.name}!")
            return

        print(f" -> {out_path.name}")

        cfg = self.plot_config

        plt.rcParams.update(cfg.depth_profile_plot_rcParams)
        fig, ax = plt.subplots(figsize=(cfg.depth_profile_plot_width, cfg.depth_profile_plot_height))

        all_clim_lows = [data["clim"][0] for data in profile_data_grouped.values() if data]
        all_clim_highs = [data["clim"][1] for data in profile_data_grouped.values() if data]
        clim = (min(all_clim_lows), max(all_clim_highs)) if all_clim_lows else (0, 1)
        cmap = plt.get_cmap("tab20")

        x_range = clim[1] - clim[0]
        expand = x_range * 0.05
        new_xlim = (clim[0] - expand, clim[1] + expand)

        ymin = max(cfg.depth_profile_surface_depth, reaction_depth - y_offset)
        ymax = reaction_depth + y_offset

        for i, (model_id, data) in enumerate(profile_data_grouped.items()):
            if not data:
                continue

            depths = np.array(data["depths"], dtype=float)
            values = np.array(data["values"], dtype=float)

            sort_idx = np.argsort(depths)
            depths = depths[sort_idx]
            values = values[sort_idx]

            unique_depths, indices = np.unique(depths, return_inverse=True)
            values_agg = np.array([values[indices == j].mean() for j in range(len(unique_depths))])

            depth_grid = np.linspace(unique_depths.min(), unique_depths.max(), depth_points)
            values_interp = np.interp(depth_grid, unique_depths, values_agg)

            dz = np.median(np.diff(depth_grid))
            window_pts = int(np.ceil(smoothing_range / dz))
            if window_pts % 2 == 0:
                window_pts += 1

            window_pts = min(window_pts, len(values_interp) if len(values_interp) % 2 == 1 else len(values_interp) - 1)
            window_pts = max(3, window_pts)

            values_smooth = savgol_filter(values_interp, window_length=window_pts, polyorder=2)

            label = model_id.replace("_", "-")
            color = cmap(i % 20)
            ax.plot(values_smooth, depth_grid / 1e3, label=label, color=color)

        ax.set_xlabel(cfg.bar_mapping.get(field_name, field_name))
        ax.set_ylabel("Depth (km)")
        ax.set_xlim(new_xlim)
        ax.set_ylim(ymin / 1e3, ymax / 1e3)
        ax.invert_yaxis()
        ax.grid(True, which="both", linewidth=0.5, color="#999999")
        ax.tick_params(axis="both", which="both", length=0)
        ax.axhline(reaction_depth / 1e3, color="black", linestyle="--", linewidth=1, zorder=1)
        ax.legend(fontsize="small", loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1)
        plt.title(f"410 {model_type.title()}")
        plt.tight_layout()
        plt.savefig(out_path, dpi=cfg.depth_profile_plot_rcParams.get("figure.dpi", 300))
        plt.close(fig)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def draw_displacement(self, displacement_data_grouped: dict[str, dict], model_type: str, out_path: Path) -> None:
        """"""
        if out_path.exists():
            print(f" -- Found plot: {out_path.name}!")
            return

        print(f" -> {out_path.name}")

        cfg = self.plot_config

        all_models = sorted(displacement_data_grouped.keys(), key=self._extract_sort_key)

        plt.rcParams.update(cfg.depth_profile_plot_rcParams)
        fig, ax = plt.subplots(figsize=(cfg.topography_plot_width, cfg.topography_plot_height))
        cmap = plt.get_cmap("tab20")

        for i, model_id in enumerate(all_models):
            data = displacement_data_grouped[model_id]
            times = data["time_myr"]
            displacements = data["displacement"]

            if len(times) > 0 and len(displacements) > 0:
                label = model_id.replace("_", "-")
                color = cmap(i % 20)
                ax.plot(times, displacements / 1e3, label=label, color=color)

        ax.set_xlabel("Time (Ma)")
        ax.set_ylabel("Displacement (km)")
        ax.grid(True, which="both", linewidth=0.5, color="#999999")
        ax.legend(fontsize="small", loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1)
        plt.title(f"410 {model_type.title()}")
        plt.tight_layout()
        plt.savefig(out_path, dpi=cfg.depth_profile_plot_rcParams.get("figure.dpi", 300))
        plt.close(fig)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def draw_width(self, width_data_grouped: dict[str, dict], model_type: str, out_path: Path) -> None:
        """"""
        if out_path.exists():
            print(f" -- Found plot: {out_path.name}!")
            return

        print(f" -> {out_path.name}")

        cfg = self.plot_config

        all_models = sorted(width_data_grouped.keys(), key=self._extract_sort_key)

        plt.rcParams.update(cfg.depth_profile_plot_rcParams)
        fig, ax = plt.subplots(figsize=(cfg.topography_plot_width, cfg.topography_plot_height))
        cmap = plt.get_cmap("tab20")

        for i, model_id in enumerate(all_models):
            data = width_data_grouped[model_id]
            times = data["time_myr"]
            widths = data["width"]

            if len(times) > 0 and len(widths) > 0:
                label = model_id.replace("_", "-")
                color = cmap(i % 20)
                ax.plot(times, widths / 1e3, label=label, color=color)

        ax.set_xlabel("Time (Ma)")
        ax.set_ylabel("Width (km)")
        ax.grid(True, which="both", linewidth=0.5, color="#999999")
        ax.legend(fontsize="small", loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1)
        plt.title(f"410 {model_type.title()}")
        plt.tight_layout()
        plt.savefig(out_path, dpi=cfg.depth_profile_plot_rcParams.get("figure.dpi", 300))
        plt.close(fig)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_ordered_pvtu_files(self) -> dict[str, tuple[list[Path], list[int]]]:
        """"""
        model_entries: list[tuple[str, tuple[list[Path], list[int]]]] = []

        for model_id, directory in self.pvtu_in_dirs.items():
            solution_dir = directory / "solution"
            if not solution_dir.is_dir():
                if self.verbosity >= 1:
                    print(f" !! Warning: solution directory not found for model: {model_id}")
                model_entries.append((model_id, ([], [])))
                continue

            try:
                files = list(solution_dir.glob("*.pvtu"))
                if not files:
                    if self.verbosity >= 1:
                        print(f" !! Warning: no .pvtu files found for model: {model_id}")
                    model_entries.append((model_id, ([], [])))
                    continue
            except Exception as e:
                if self.verbosity >= 1:
                    print(f" !! Warning: error processing directory {solution_dir}:\n    {e}")
                model_entries.append((model_id, ([], [])))
                continue

            parsed = []
            for file in files:
                match = re.search(r"solution-(\d+)\.pvtu", file.name)
                if match:
                    timestep = int(match.group(1))
                else:
                    timestep = -1
                    if self.verbosity >= 1:
                        print(f" !! Warning: could not parse timestep from filename:\n    {file.name}")
                parsed.append((timestep, file))

            # Sort by timestep
            parsed.sort(key=lambda x: x[0])
            sorted_timesteps = [t for t, _ in parsed]
            sorted_files = [f for _, f in parsed]

            # Sort by model id
            model_entries.append((model_id, (sorted_files, sorted_timesteps)))

        model_entries.sort(key=lambda x: self._extract_sort_key(x[0]))

        # Convert back to dict to preserve sorted order
        ordered_pvtu_files: dict[str, tuple[list[Path], list[int]]] = dict(model_entries)

        return ordered_pvtu_files

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _extract_z_b_values(self, model_id_str) -> tuple[float, float]:
        """"""
        try:
            parts = model_id_str.split("_")

            Z_value_str = parts[1][1:]
            B_value_str = parts[2][1:]

            Z_factor = float(Z_value_str)
            B_factor = float(B_value_str)

            return Z_factor, B_factor
        except (IndexError, ValueError) as e:
            print(f" !! Warning: could not extract Z/B values from model_id '{model_id_str}':\n    {e}")
            return np.nan, np.nan

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _extract_sort_key(self, model_id: str):
        """"""
        match = re.match(r"(.*?)(\d+)$", model_id)
        if match:
            prefix, number = match.groups()
            return (prefix, int(number))
        else:
            return (model_id, float("inf"))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _prepare_mesh(self, mesh: pv.DataSet, field_name: str) -> None:
        """Ensures necessary fields are on the mesh, calculating them if possible."""
        if re.fullmatch(r"stress_[xyz]{2}", field_name) and "stress" in mesh.point_data:
            if field_name not in mesh.point_data:
                idx = self.plot_config.stress_components.get(field_name)
                if idx is not None:
                    mesh[field_name] = mesh["stress"][:, idx]

        elif re.fullmatch(r"shear_stress_[xyz]{2}", field_name) and "shear_stress" in mesh.point_data:
            if field_name not in mesh.point_data:
                idx = self.plot_config.shear_stress_components.get(field_name)
                if idx is not None:
                    mesh[field_name] = mesh["shear_stress"][:, idx]

        elif field_name == "differential_stress":
            if "principal_stress_1" in mesh.point_data and "principal_stress_2" in mesh.point_data:
                if field_name not in mesh.point_data:
                    mesh[field_name] = mesh["principal_stress_1"] - mesh["principal_stress_2"]
            else:
                if self.verbosity >= 1:
                    print(" !! Warning: cannot calculate 'differential_stress': missing principal stress components.")

        elif field_name == "nonadiabatic_density":
            if "density" in mesh.point_data and "adiabatic_density" in mesh.point_data:
                if field_name not in mesh.point_data:
                    mesh[field_name] = mesh["density"] - mesh["adiabatic_density"]
            else:
                if self.verbosity >= 1:
                    print(" !! Warning: cannot calculate 'nonadiabatic_density': missing density components.")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _check_mesh(self, mesh: pv.UnstructuredGrid, field_name: str) -> bool:
        """"""
        if field_name not in mesh.point_data and field_name not in self.plot_config.get_extra_fields():
            return False

        if field_name not in mesh.point_data:
            return False

        return True

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_mesh_time_myr(self, mesh: pv.UnstructuredGrid, model_id: str, tstep: int) -> float:
        """"""
        if mesh.field_data is not None and "TIME" in mesh.field_data:
            time_myr = mesh.field_data["TIME"][0] / 1e6
        else:
            time_myr = 0.0
            print(f"'TIME' field not found in mesh.field_data for timestep {tstep}, model {model_id}.")

        return time_myr

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _find_vertical_profile_for_410_structure(
        self, mesh: pv.UnstructuredGrid, field_name: str, reaction_depth: float, Z_factor: float, B_factor: float, width_threshold: float = 15e3
    ) -> tuple[float, float, float]:
        """"""
        cfg = self.plot_config

        surface_depth = cfg.depth_profile_surface_depth

        if (field_name == "velocity" or cfg.velocity_mapping.get(field_name, False)) and "velocity" in mesh.point_data:
            mesh["velocity"] = np.abs(mesh["velocity"][:, 1])

        x_offset = (
            50e3
            if (B_factor <= 4 and Z_factor < 2.0e2) or (B_factor > 4 and B_factor <= 6 and Z_factor < 3.7e1) or (B_factor > 6 and Z_factor < 1.6e1)
            else 450e3
        )

        x_coords = mesh.points[:, 0]
        x_min, x_max = x_coords.min(), x_coords.max()
        x_center = 0.5 * (x_min + x_max)
        x_start = max(x_center - x_offset, x_min)
        x_end = min(x_center + x_offset, x_max)
        all_unique_x = np.unique(x_coords)
        x_positions = all_unique_x[(all_unique_x >= x_start) & (all_unique_x <= x_end)]

        results = []

        for x in x_positions:
            mask = x_coords == x
            if not np.any(mask):
                continue

            depths = mesh.point_data["depth"][mask] + surface_depth
            values = mesh.point_data[field_name][mask]

            sort_idx = np.argsort(depths)
            depths = depths[sort_idx]
            values = values[sort_idx]

            displacement, width = self._measure_phase_transition_displacement_and_width(depths, values, reaction_depth)
            if np.isnan(displacement) or np.isnan(width):
                continue

            results.append((x, displacement, width))

        if not results:
            return np.nan, np.nan, np.nan

        results = np.array(results)
        x_vals, displacements, widths = results.T

        measure_displacement = (
            (B_factor >= 6 and B_factor < 10 and Z_factor > 4.3e5)
            or (B_factor >= 10 and Z_factor > 6.0e3)
            or (B_factor < 6 and Z_factor < 8.7e1)
            or (B_factor <= 10 and Z_factor <= 1.6e1)
        )

        if np.any(widths > width_threshold) and not measure_displacement:
            best_idx = np.nanargmax(widths)
        else:
            best_idx = np.nanargmax(np.abs(displacements))

        best_x = x_vals[best_idx]
        best_displacement = displacements[best_idx]
        best_width = widths[best_idx]

        return best_displacement, best_width, best_x

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _measure_phase_transition_displacement_and_width(self, depths: np.ndarray, values: np.ndarray, reaction_depth: float) -> tuple[float, float]:
        """"""

        def get_first_crossing_index(target: float, start_index: int = 0) -> int:
            """"""
            if start_index >= len(values) - 1:
                return -1

            v = values[start_index:]
            diffs = v - target

            sign_change_indices = np.where(np.diff(np.sign(diffs)) != 0)[0]

            if not sign_change_indices.any():
                return -1

            return sign_change_indices[0] + start_index

        def interpolate_at_index(target: float, i: int) -> float:
            """"""
            if i < 0 or i >= len(values) - 1:
                return np.nan

            x0, x1 = values[i], values[i + 1]
            z0, z1 = depths[i], depths[i + 1]

            if x1 == x0:
                return float((z0 + z1) / 2)

            weight = (target - x0) / (x1 - x0)
            return float(z0 + weight * (z1 - z0))

        if len(depths) < 2:
            return np.nan, np.nan

        idx_10 = get_first_crossing_index(0.1, 0)
        idx_90 = -1

        if idx_10 != -1:
            if values[idx_10 + 1] >= 1.0:
                idx_90 = idx_10
            else:
                idx_90 = get_first_crossing_index(0.9, idx_10 + 1)

        if idx_10 == -1 or idx_90 == -1:
            return np.nan, np.nan

        depth_at_X10 = interpolate_at_index(0.1, idx_10)
        depth_at_X90 = interpolate_at_index(0.9, idx_90)

        if np.isnan(depth_at_X10) or np.isnan(depth_at_X90):
            return np.nan, np.nan

        displacement = depth_at_X90 - reaction_depth
        width = depth_at_X90 - depth_at_X10

        return displacement, width

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_depth_profile_data_at_x_pos(self, mesh: pv.UnstructuredGrid, field_name: str, x_pos: float) -> tuple[np.ndarray, np.ndarray]:
        """"""
        cfg = self.plot_config

        surface_depth = cfg.depth_profile_surface_depth

        x_coords = mesh.points[:, 0]
        mask = x_coords == x_pos

        if not np.any(mask):
            return np.empty(0), np.empty(0)

        if (field_name == "velocity" or cfg.velocity_mapping.get(field_name, False)) and "velocity" in mesh.point_data:
            mesh["velocity"] = np.abs(mesh["velocity"][:, 1])

        depths = mesh.point_data["depth"][mask] + surface_depth
        values = mesh.point_data[field_name][mask] * cfg.scale_mapping.get(field_name, 1.0)

        sort_idx = np.argsort(depths)
        return depths[sort_idx], values[sort_idx]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _is_diverging_cmap(self, cmap_name: str) -> bool:
        diverging_base_names = [
            "PiYG",
            "PRGn",
            "BrBG",
            "PuOr",
            "RdBu",
            "RdGy",
            "RdYlBu",
            "RdYlGn",
            "Spectral",
            "coolwarm",
            "bwr",
            "seismic",
        ]
        diverging_names = diverging_base_names + [name + "_r" for name in diverging_base_names]

        return cmap_name in diverging_names

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _is_sequential_cmap(self, cmap_name: str) -> bool:
        sequential_base_names = [
            "viridis",
            "plasma",
            "inferno",
            "magma",
            "cividis",
            "gist_heat",
            "pink",
            "bone",
            "afmhot",
            "copper",
            "Greys",
            "Purples",
            "Blues",
            "Greens",
            "Oranges",
            "Reds",
            "YlOrBr",
            "YlOrRd",
            "OrRd",
            "PuRd",
            "RdPu",
            "BuPu",
            "GnBu",
            "PuBu",
            "YlGnBu",
            "PuBuGn",
            "BuGn",
            "YlGn",
        ]

        sequential_names = sequential_base_names + [name + "_r" for name in sequential_base_names]

        return cmap_name in sequential_names

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _modify_diverging_cmap(self, cmap_name: str, central_color_hex: str, n_colors: int = 11) -> mcolors.Colormap | mcolors.ListedColormap:
        if not self._is_diverging_cmap(cmap_name):
            print(f" !! Warning: '{cmap_name}' is not a recognized diverging colormap!")
            return plt.get_cmap(cmap_name, n_colors)

        original_cmap = plt.get_cmap(cmap_name)
        cmap_colors = original_cmap(np.linspace(0, 1, n_colors))
        center_idx = n_colors // 2
        cmap_colors[center_idx] = mcolors.to_rgba(central_color_hex)

        return mcolors.ListedColormap(cmap_colors)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _modify_sequential_cmap(self, cmap_name: str, end_color_hex: str, n_colors: int = 11) -> mcolors.Colormap | mcolors.ListedColormap:
        if not self._is_sequential_cmap(cmap_name):
            print(f" !! Warning: '{cmap_name}' is not a recognized sequential colormap!")
            return plt.get_cmap(cmap_name, n_colors)

        original_cmap = plt.get_cmap(cmap_name)

        cmap_colors = original_cmap(np.linspace(0, 1, n_colors))

        change_first_color_cmaps = [
            "gist_heat_r",
            "pink_r",
            "bone_r",
            "afmhot_r",
            "copper_r",
            "Purples",
            "BuPu",
            "Reds",
            "Oranges",
            "YlOrBr",
            "GnBu",
        ]

        if cmap_name in change_first_color_cmaps:
            cmap_colors[0] = mcolors.to_rgba(end_color_hex)
        else:
            cmap_colors[-1] = mcolors.to_rgba(end_color_hex)

        return mcolors.ListedColormap(cmap_colors)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _add_X_field_contours(
        self, plotter: pv.Plotter, mesh: pv.UnstructuredGrid, color: str = "black", alpha: float = 1.0, line_width: int = 5
    ) -> None:
        """"""

        if "X_field" not in mesh.point_data:
            return

        iso_values = [0.1, 0.9]

        try:
            contour = mesh.contour(isosurfaces=iso_values, scalars="X_field")
        except Exception as e:
            print(f" !! Warning: failed generating X_field contours:\n    {e}")
            return

        plotter.add_mesh(contour, color=color, opacity=alpha, line_width=line_width, render_lines_as_tubes=True)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _add_vertical_contour(
        self,
        plotter: pv.Plotter,
        mesh: pv.UnstructuredGrid,
        x_pos: float,
        color: str = "black",
        alpha: float = 1.0,
        line_width: int = 5,
        tolerance_km: float = 2,
        verbosity: int = 0,
    ) -> None:
        """"""

        x_coords = mesh.points[:, 0]

        mask = np.isclose(x_coords, x_pos, atol=tolerance_km * 1e3)
        contour_points = mesh.points[mask]

        if contour_points.shape[0] > 1:
            y_coords = contour_points[:, 1]
            z_coords = contour_points[:, 2] if contour_points.shape[1] > 2 else np.array([0])

            if np.ptp(y_coords) >= np.ptp(z_coords):
                sort_idx = np.argsort(y_coords)
            else:
                sort_idx = np.argsort(z_coords)

            ordered_contour_points = contour_points[sort_idx]

            polyline = pv.lines_from_points(ordered_contour_points, close=False)

            plotter.add_mesh(polyline, color=color, line_width=line_width, opacity=alpha, render_lines_as_tubes=True)
        elif verbosity >= 1:
            print(f" !! Warning: Fewer than 2 points found at x={x_pos / 1e3:.1f} km (tolerance={tolerance_km} m).\n -- Skipping X contour")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _add_depth_contour(
        self,
        plotter: pv.Plotter,
        mesh: pv.UnstructuredGrid,
        depth: float,
        color: str = "black",
        alpha: float = 1.0,
        line_width: int = 5,
        tolerance_km: float = 2,
        verbosity: int = 1,
    ) -> None:
        """"""
        if "depth" not in mesh.point_data:
            print("'depth' field not found in mesh.\n -- Skipping depth contour")
            return

        depth_values_m = mesh.point_data["depth"]
        mask = np.isclose(depth_values_m, depth, atol=tolerance_km * 1e3)
        contour_points = mesh.points[mask]

        if contour_points.shape[0] > 1:
            x_min, x_max = np.min(contour_points[:, 0]), np.max(contour_points[:, 0])
            y_val = np.mean(contour_points[:, 1])

            line_points = np.array([[x_min, y_val, 0.0], [x_max, y_val, 0.0]])
            polyline = pv.lines_from_points(line_points)
            plotter.add_mesh(polyline, color=color, line_width=line_width, opacity=alpha, render_lines_as_tubes=True)
        elif verbosity >= 1:
            print(f" !! Warning: Fewer than 2 points found at depth {depth / 1e3:.1f} km (tolerance={tolerance_km} m).\n -- Skipping X contour")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _configure_cmap(
        self, mesh: pv.UnstructuredGrid, field_name: str, expansion: float = 0.0
    ) -> tuple[mcolors.Colormap | mcolors.ListedColormap, tuple[float, float]]:
        """"""
        plot_config = self.plot_config

        cmap_choice = plot_config.cmap_mapping.get(field_name, "viridis")
        current_scalars = mesh.point_data[field_name] * plot_config.scale_mapping.get(field_name, 1.0)

        if field_name in ["viscosity", "strain_rate", "arrhenius"]:
            current_scalars = np.log10(np.abs(np.maximum(current_scalars, 1e-30)))
        elif field_name == "reaction_rate_C0":
            current_scalars = np.sign(current_scalars) * np.log10(1.0 + np.abs(current_scalars))
            cmap_choice = plot_config.cmap_mapping.get(field_name, "seismic")

        mesh[f"{field_name}_viz"] = current_scalars

        clim_config = plot_config.clim_mapping.get(field_name, None)
        if clim_config in (None, "auto"):
            finite_data = current_scalars[np.isfinite(current_scalars)]
            if finite_data.size > 0:
                clim_min = np.min(finite_data)
                clim_max = np.max(finite_data)
                if clim_min == clim_max:
                    clim_actual = (clim_min, clim_max)
                else:
                    if self._is_diverging_cmap(cmap_choice):
                        clim_max_abs = max(abs(clim_min), abs(clim_max))
                        clim_actual = (-clim_max_abs, clim_max_abs)
                    else:
                        clim_actual = (clim_min, clim_max)
            else:
                clim_actual = (0, 1)
            if self.verbosity >= 1:
                print(f" !! Warning: Using auto CLIM for '{field_name}': {clim_actual}")
        else:
            if isinstance(clim_config, (tuple, list)) and len(clim_config) == 2:
                clim_actual = tuple(map(float, clim_config))
            else:
                clim_actual = (0, 1)

        if expansion != 0.0:
            try:
                cmin, cmax = float(clim_actual[0]), float(clim_actual[1])
                crange = cmax - cmin
                clim_actual = (cmin - expansion * crange, cmax + expansion * crange)
            except (TypeError, ValueError) as e:
                if self.verbosity >= 1:
                    print(f" !! Warning: expansion skipped: CLIM values not numeric: {clim_actual}:\n    {e}")

        modified_seq_target = ["stress_second_invariant", "X_field"]

        if self._is_diverging_cmap(cmap_choice):
            central_color = "#FFFFFF" if cmap_choice in ["RdGy", "RdGy_r"] else "#E5E5E5"
            final_cmap = self._modify_diverging_cmap(cmap_choice, central_color, plot_config.n_colors)
        elif self._is_sequential_cmap(cmap_choice) and field_name in modified_seq_target:
            final_cmap = self._modify_sequential_cmap(cmap_choice, "#E5E5E5", plot_config.n_colors)
        else:
            final_cmap = plt.get_cmap(cmap_choice, plot_config.n_colors)

        return final_cmap, clim_actual

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _compute_camera_settings(self, mesh, full_view: bool) -> list[tuple[float, ...]]:
        """"""
        bounds = mesh.bounds

        center_pt = [
            (bounds[0] + bounds[1]) / 2,
            (bounds[2] + bounds[3]) / 2,
            (bounds[4] + bounds[5]) / 2,
        ]
        max_bound_size = max(bounds[1] - bounds[0], bounds[3] - bounds[2], 1e-6)

        if not self.plot_config.camera_center_zoom:
            if full_view:
                camera_position = [
                    (0, 0, bounds[1] * 3.8),
                    (0, 0, 0),
                    (0, 1, 0),
                ]
            else:
                camera_position = [
                    (bounds[1] * 0, bounds[1] + 0.67, bounds[1] * 1.3),
                    (bounds[1] * 0, bounds[1] + 0.67, 0),
                    (0, 1, 0),
                ]
        else:
            cam_y_shift = max_bound_size * self.plot_config.camera_y_shift_factor
            cam_dist = max_bound_size * self.plot_config.camera_zoom_factor

            camera_position = [
                (
                    center_pt[0],
                    center_pt[1],
                    cam_dist if cam_dist > 0 else 1000,
                ),
                (center_pt[0], center_pt[1] + cam_y_shift, center_pt[2]),
                (0, 1, 0),
            ]

        return camera_position

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _add_scale_bar(self, mesh, plotter) -> None:
        """"""
        # Add scale bar
        bounds = mesh.bounds
        x_range = bounds[1] - bounds[0]
        scale_length = x_range * self.plot_config.scale_bar_length_fraction
        scale_km = scale_length / 1000

        start_point = [
            bounds[0] + x_range * self.plot_config.scale_bar_position[0],
            bounds[2] + (bounds[3] - bounds[2]) * self.plot_config.scale_bar_position[1],
            0,
        ]
        end_point = [start_point[0] + scale_length, start_point[1], 0]

        line = pv.Line(start_point, end_point)
        plotter.add_mesh(line, color=self.plot_config.scale_bar_color, line_width=self.plot_config.scale_bar_thickness)

        text_pos = [start_point[0], start_point[1] + 0.01 * x_range, 0]
        shift_factor = self.plot_config.scale_bar_label_shift_factor
        text_pos[0] += shift_factor * x_range

        plotter.add_point_labels(
            [text_pos],
            [f"{scale_km:.0f} km"],
            font_size=self.plot_config.scale_bar_label_font_size,
            point_color=self.plot_config.scale_bar_color,
            text_color=self.plot_config.scale_bar_color,
            point_size=0,
            shape=None,
        )
