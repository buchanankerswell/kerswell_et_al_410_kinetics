#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
import gc
import glob
import re
import warnings
from argparse import ArgumentParser, Namespace
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from PIL import Image
from scipy.spatial import ConvexHull

__import__("vtk")


#######################################################
## .1. PlottingConfig                            !!! ##
#######################################################
@dataclass
class PlottingConfig:
    """Holds all configuration for PyVista visualizations."""

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
            "Vp_Vs_ratio",
        ]
    )

    n_colors: int = 11
    show_edges: bool = False
    default_fig_dir: Path | str = "figs/simulation"
    plotter_window_size: list[int] = field(default_factory=lambda: [1920, 1920])
    plotter_background: str = "#FFFFFF"
    plotter_lighting: str = "none"
    edge_opacity: float = 0.3
    title_font_size: int = 60
    title_position: tuple[float, float] = field(default_factory=lambda: (0.39, 0.40))
    cbar_vertical: bool = False
    cbar_title_font_size: int = 115
    cbar_label_font_size: int = 80
    cbar_width: float = 0.4
    cbar_n_labels: int = 3
    cbar_position: list[float] = field(default_factory=lambda: [0.3, 0.5])
    rotation_init: float = 0
    rotation_coeff: float = 0
    camera_center_zoom: bool = False
    camera_full_view: bool = True
    camera_y_shift_factor: float = -0.045
    camera_zoom_factor: float = 1.9
    screenshot_dpi: tuple[int, int] = field(default_factory=lambda: (330, 330))
    depth_contour_depths_km: list[int] = field(
        default_factory=lambda: [0, 125, 410, 660, 2890]
    )
    depth_contour_line_widths: list[int] = field(
        default_factory=lambda: [6, 3, 3, 3, 6]
    )
    depth_contour_tolerance_km: float = 15
    velocity_glyph_factor: float = 7e5
    stress_glyph_factor: float = 1e-3
    glyph_line_width: float = 1
    scale_bar_enabled: bool = True
    scale_bar_color: str = "black"
    scale_bar_thickness: float = 25
    scale_bar_length_fraction: float = 0.25
    scale_bar_position: tuple = (0.05, 0.08)
    scale_bar_font_size_factor: float = 1.0

    centerline_depth_plots: bool = False
    centerline_plot_fig_width: float = 4.5
    centerline_plot_fig_height: float = 6.5
    centerline_plot_rcParams: dict[str, Any] = field(
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
                "v_r": "velocity-rad",
                "v_phi": "velocity-phi",
                "principal_stress_1": "sigma-one",
                "principal_stress_2": "sigma-two",
                "differential_stress": "sigma-delta",
                "seismic_Vp": "vp",
                "seismic_Vs": "vs",
                "Vp_anomaly": "vp-anomaly",
                "Vs_anomaly": "vs-anomaly",
                "Vp_Vs_ratio": "vp-vs-ratio",
                "driving_force": "driving-force",
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
                "v_r": "BrBG",
                "v_phi": "BrBG",
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
                "seismic_Vp": "bone_r",
                "seismic_Vs": "bone_r",
                "Vp_anomaly": "seismic",
                "Vs_anomaly": "seismic",
                "Vp_Vs_ratio": "seismic",
                "driving_force": "RdGy",
                "reaction_rate_C0": "bone_r",
            }
        if not self.clim_mapping:
            self.clim_mapping = {
                "T": (273, 3773),
                "adiabatic_temperature": (273, 3773),
                "nonadiabatic_temperature": (-1405, 1405),
                "p": (-0.5, 136.5),
                "adiabatic_pressure": (-0.5, 136.5),
                "nonadiabatic_pressure": (-1.85, 1.85),
                "X_field": (0.0, 1.0),
                "density": (3.0, 5.8),
                "adiabatic_density": (3.0, 5.8),
                "nonadiabatic_density": (-0.29, 0.29),
                "adiabatic_density_derivative": (0, 12),
                "specific_heat": (1.15, 1.33),
                "thermal_expansivity": (1.10, 4.32),
                "compressibility": (1.50e-12, 5.35e-11),
                "viscosity": (18.5, 25),
                "velocity": (-20, 20),
                "v_r": (-20, 20),
                "v_phi": (-20, 20),
                "stress_xx": (0, 136.5),
                "stress_yy": (0, 136.5),
                "shear_stress_xx": (-1.58, 1.58),
                "shear_stress_xy": (-1.58, 1.58),
                "shear_stress_yy": (-1.58, 1.58),
                "principal_stress_1": (-1.0, 1.0),
                "principal_stress_2": (-1.0, 1.0),
                "differential_stress": (0, 2.0),
                "stress_second_invariant": (0, 1.1),
                "strain_rate": (-19.1, -12.0),
                "seismic_Vp": (0, 10),
                "seismic_Vs": (0, 10),
                "Vp_anomaly": (-1, 1),
                "Vs_anomaly": (-1, 1),
                "Vp_Vs_ratio": (-4.5, 4.5),
                "driving_force": (-15, 15),
                "reaction_rate_C0": (0, 10e-14),
            }
        if not self.bar_mapping:
            self.bar_mapping = {
                "T": "$T$ (K)",
                "adiabatic_temperature": "$\\bar{T}$ (K)",
                "nonadiabatic_temperature": "$\\hat{T}$ (K)",
                "p": "$P$ (GPa)",
                "adiabatic_pressure": "$\\bar{P}$ (GPa)",
                "nonadiabatic_pressure": "$\\hat{P}$ (GPa)",
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
                "v_r": "$\\vec{u}_r$ (cm/yr)",
                "v_phi": "$\\vec{u}_{\\phi}$ (cm/yr)",
                "stress_xx": "$\\sigma_{xx}$ (GPa)",
                "stress_yy": "$\\sigma_{yy}$ (GPa)",
                "shear_stress_xx": "$\\sigma^{\\prime}_{xx}$ (GPa)",
                "shear_stress_xy": "$\\sigma^{\\prime}_{xy}$ (GPa)",
                "shear_stress_yy": "$\\sigma^{\\prime}_{yy}$ (GPa)",
                "principal_stress_1": "$\\sigma^{\\prime}_1$ (GPa)",
                "principal_stress_2": "$\\sigma^{\\prime}_2$ (GPa)",
                "differential_stress": "$\\sigma^{\\prime}_1$-$\\sigma^{\\prime}_2$ (GPa)",
                "stress_second_invariant": "$\\sigma_{II}$ (GPa)",
                "strain_rate": "Log $\\dot{\\epsilon}_{II}$ (1/s)",
                "seismic_Vp": "$V_p$ (km/s)",
                "seismic_Vs": "$V_s$ (km/s)",
                "Vp_anomaly": "$\\hat{V}_p$ (%)",
                "Vs_anomaly": "$\\hat{V}_s$ (%)",
                "Vp_Vs_ratio": "$V_p$/$V_s$",
                "driving_force": "$\\Delta G (kJ/mol)$",
                "reaction_rate_C0": "$dX/dt (1/Ma)$",
            }
        if not self.fmt_mapping:
            self.fmt_mapping = {
                "T": "%.0f",
                "adiabatic_temperature": "%.0f",
                "nonadiabatic_temperature": "%.0f",
                "p": "%.0f",
                "adiabatic_pressure": "%.0f",
                "nonadiabatic_pressure": "%.1f",
                "X_field": "%.1f",
                "density": "%.1f",
                "adiabatic_density": "%.1f",
                "nonadiabatic_density": "%.2f",
                "adiabatic_density_derivative": "%.0f",
                "specific_heat": "%.2f",
                "thermal_expansivity": "%.1f",
                "compressibility": "%.1e",
                "viscosity": "%.0f",
                "velocity": "%.0f",
                "v_r": "%.0f",
                "v_phi": "%.0f",
                "stress_xx": "%.0f",
                "stress_yy": "%.0f",
                "shear_stress_xx": "%.1f",
                "shear_stress_xy": "%.1f",
                "shear_stress_yx": "%.1f",
                "shear_stress_yy": "%.1f",
                "principal_stress_2": "%.1f",
                "principal_stress_1": "%.1f",
                "differential_stress": "%.1f",
                "stress_second_invariant": "%.1f",
                "strain_rate": "%.0f",
                "seismic_Vp": "%.1f",
                "seismic_Vs": "%.1f",
                "Vp_anomaly": "%.0f",
                "Vs_anomaly": "%.0f",
                "Vp_Vs_ratio": "%.0f",
                "driving_force": "%.0f",
                "reaction_rate_C0": "%.1f",
            }
        if not self.scale_mapping:
            self.scale_mapping = {
                "T": 1,
                "adiabatic_temperature": 1,
                "nonadiabatic_temperature": 1,
                "p": 1e-9,
                "adiabatic_pressure": 1e-9,
                "nonadiabatic_pressure": 1e-9,
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
                "v_r": 1e2,
                "v_phi": 1e2,
                "stress_xx": 1e-9,
                "stress_yy": 1e-9,
                "shear_stress_xx": 1e-9,
                "shear_stress_xy": 1e-9,
                "shear_stress_yy": 1e-9,
                "principal_stress_1": 1e-9,
                "principal_stress_2": 1e-9,
                "differential_stress": 1e-9,
                "stress_second_invariant": 1e-9,
                "strain_rate": 1,
                "seismic_Vp": 1e-3,
                "seismic_Vs": 1e-3,
                "Vp_anomaly": 1,
                "Vs_anomaly": 1,
                "Vp_Vs_ratio": 1,
                "driving_force": 1e-3,
                "reaction_rate_C0": 3.154e13,
            }
        if not self.velocity_mapping:
            self.velocity_mapping = {
                "T": False,
                "adiabatic_temperature": False,
                "nonadiabatic_temperature": False,
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
                "v_r": False,
                "v_phi": False,
                "principal_stress_1": False,
                "principal_stress_2": False,
                "differential_stress": False,
                "seismic_Vp": False,
                "seismic_Vs": False,
                "Vp_anomaly": False,
                "Vs_anomaly": False,
                "Vp_Vs_ratio": False,
                "driving_force": False,
                "reaction_rate_C0": False,
            }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def get_extra_fields(self) -> list[str]:
        """Returns a list of fields that can be derived from other mesh fields."""
        return (
            list(self.stress_components.keys())
            + list(self.shear_stress_components.keys())
            + self.add_mesh_fields
        )


#######################################################
## .2. PyVistaModelVisualizer                    !!! ##
#######################################################
@dataclass
class PyVistaModelVisualizer:
    plot_config: "PlottingConfig"
    pvtu_in_dir: dict[str, Path]
    tsteps: list[int] | None
    verbosity: int = 0

    def __post_init__(self):
        if self.tsteps is None:
            self.tsteps = []

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_pvtu_files(self) -> dict[str, list[str]]:
        """"""
        pvtu_files: dict[str, list[str]] = {}
        for model_id, directory in self.pvtu_in_dir.items():
            solution_dir = directory / "solution"
            if not solution_dir.is_dir():
                if self.verbosity >= 1:
                    print(" !! Warning: solution directory not found!")
                pvtu_files[model_id] = []
                continue

            try:
                files = glob.glob(str(solution_dir / "*.pvtu"))
                if not files:
                    if self.verbosity >= 1:
                        print(" !! Warning: no .pvtu files found!")
                    pvtu_files[model_id] = []
                    continue

                def extract_step_number(filename: Path | str) -> int:
                    if isinstance(filename, str):
                        filename = Path(filename)
                    match = re.search(r"solution-(\d+)\.pvtu", filename.name)
                    if match:
                        return int(match.group(1))
                    else:
                        if self.verbosity >= 1:
                            print(
                                f" !! Warning: filename does not match expected pattern: {filename}"
                            )
                        return -1

                pvtu_files[model_id] = sorted(files, key=extract_step_number)

            except Exception as e:
                if self.verbosity >= 1:
                    print(
                        f" !! Warning: error processing directory {solution_dir}:\n    {e}"
                    )
                pvtu_files[model_id] = []

        return pvtu_files

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _prepare_mesh_fields(self, mesh: pv.DataSet, field_name: str) -> pv.DataSet:
        """Ensures necessary fields are on the mesh, calculating them if possible."""
        if re.fullmatch(r"stress_[xyz]{2}", field_name) and "stress" in mesh.point_data:
            if field_name not in mesh.point_data:
                idx = self.plot_config.stress_components.get(field_name)
                if idx is not None:
                    mesh[field_name] = mesh["stress"][:, idx]

        elif (
            re.fullmatch(r"shear_stress_[xyz]{2}", field_name)
            and "shear_stress" in mesh.point_data
        ):
            if field_name not in mesh.point_data:
                idx = self.plot_config.shear_stress_components.get(field_name)
                if idx is not None:
                    mesh[field_name] = mesh["shear_stress"][:, idx]

        elif field_name == "differential_stress":
            if (
                "principal_stress_1" in mesh.point_data
                and "principal_stress_2" in mesh.point_data
            ):
                if field_name not in mesh.point_data:
                    mesh[field_name] = (
                        mesh["principal_stress_1"] - mesh["principal_stress_2"]
                    )
            else:
                if self.verbosity >= 1:
                    print(
                        " !! Warning: cannot calculate 'differential_stress': missing principal stress components."
                    )

        elif field_name == "nonadiabatic_density":
            if "density" in mesh.point_data and "adiabatic_density" in mesh.point_data:
                if field_name not in mesh.point_data:
                    mesh[field_name] = mesh["density"] - mesh["adiabatic_density"]
            else:
                if self.verbosity >= 1:
                    print(
                        " !! Warning: cannot calculate 'nonadiabatic_density': missing density components."
                    )

        elif field_name == "Vp_Vs_ratio":
            if "Vp_anomaly" in mesh.point_data and "Vs_anomaly" in mesh.point_data:
                if field_name not in mesh.point_data:
                    vp_anom = mesh["Vp_anomaly"]
                    vs_anom = mesh["Vs_anomaly"]
                    mesh[field_name] = np.divide(
                        vp_anom, vs_anom, out=np.zeros_like(vp_anom), where=vs_anom != 0
                    )
            else:
                if self.verbosity >= 1:
                    print(
                        " !! Warning: cannot calculate 'Vp_Vs_ratio': missing Vp/Vs anomaly data."
                    )

        return mesh

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _rotate_mesh(self, mesh: pv.DataSet, tstep: int) -> pv.DataSet:
        """"""
        plot_config = self.plot_config

        if plot_config.rotation_coeff != 0 and plot_config.rotation_init != 0:
            rotation_angle = (
                int(tstep) * plot_config.rotation_coeff
            ) + plot_config.rotation_init
            mesh.rotate_z(rotation_angle, inplace=True)

        return mesh

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _check_mesh_field(
        self, mesh: pv.DataSet, field_name: str, model_id: str, tstep: int
    ) -> bool:
        """"""
        plot_config = self.plot_config

        if (
            field_name not in mesh.point_data
            and field_name not in plot_config.get_extra_fields()
        ):
            print(
                f" -- Skipping plot: field '{field_name}' not available on mesh for timestep {tstep}, model {model_id}!"
            )
            return False

        if field_name not in mesh.point_data:
            print(
                f" -- Skipping plot: field '{field_name}' could not be prepared/calculated for mesh!"
            )
            return False

        return True

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_mesh_time_myr(self, mesh: pv.DataSet, model_id: str, tstep: int) -> float:
        """"""
        if mesh.field_data is not None and "TIME" in mesh.field_data:
            time_myr = mesh.field_data["TIME"][0] / 1e6
        else:
            time_myr = 0.0
            print(
                f"'TIME' field not found in mesh.field_data for timestep {tstep}, model {model_id}."
            )

        return time_myr

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
        diverging_names = diverging_base_names + [
            name + "_r" for name in diverging_base_names
        ]

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

        sequential_names = sequential_base_names + [
            name + "_r" for name in sequential_base_names
        ]

        return cmap_name in sequential_names

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _modify_diverging_cmap(
        self, cmap_name: str, central_color_hex: str, n_colors: int = 11
    ) -> mcolors.Colormap | mcolors.ListedColormap:
        if not self._is_diverging_cmap(cmap_name):
            print(
                f"'{cmap_name}' is not a recognized diverging colormap. Returning original."
            )
            return plt.get_cmap(cmap_name, n_colors)

        original_cmap = plt.get_cmap(cmap_name)
        cmap_colors = original_cmap(np.linspace(0, 1, n_colors))
        center_idx = n_colors // 2
        cmap_colors[center_idx] = mcolors.to_rgba(central_color_hex)

        return mcolors.ListedColormap(cmap_colors)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _modify_sequential_cmap(
        self, cmap_name: str, end_color_hex: str, n_colors: int = 11
    ) -> mcolors.Colormap | mcolors.ListedColormap:
        if not self._is_sequential_cmap(cmap_name):
            print(
                f"'{cmap_name}' is not a recognized sequential colormap. Returning original."
            )
            return plt.get_cmap(cmap_name, n_colors)

        original_cmap = plt.get_cmap(cmap_name)

        cmap_colors = original_cmap(np.linspace(0, 1, n_colors))

        change_first_color_cmaps = [
            "gist_heat_r",
            "pink_r",
            "bone_r",
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
    def _add_depth_contour(
        self,
        plotter: pv.Plotter,
        mesh: pv.DataSet,
        depth_km: float,
        color: str = "black",
        alpha: float = 1.0,
        line_width: int = 3,
        tolerance_km: float = 15,
        verbosity: int = 0,
    ) -> None:
        """Adds a depth contour line to the PyVista plotter."""
        if "depth" not in mesh.point_data:
            print("'depth' field not found in mesh.\n -- Skipping depth contour")
            return

        depth_values_m = mesh.point_data["depth"]
        mask = np.isclose(depth_values_m, depth_km * 1e3, atol=tolerance_km * 1e3)
        contour_points = mesh.points[mask]

        if contour_points.shape[0] > 2:
            xy = contour_points[:, :2]

            # Check if points are nearly colinear or have no spread in 2D
            x_range = np.ptp(xy[:, 0])
            y_range = np.ptp(xy[:, 1])

            if x_range < 1e-3 or y_range < 1e-3:
                if verbosity >= 1:
                    print(
                        f" !! Warning: points at {depth_km} km are effectively 1D!\n"
                        " -- Skipping depth contour"
                    )
                return

            try:
                hull = ConvexHull(xy)
                ordered_contour_points = contour_points[hull.vertices]
                polyline = pv.lines_from_points(ordered_contour_points, close=True)
                plotter.add_mesh(
                    polyline, color=color, line_width=line_width, opacity=alpha
                )
            except Exception as e:
                print(
                    f" !! Error: could not generate convex hull for depth contour at {depth_km}km:\n    {e}"
                )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _configure_cmap(
        self, mesh: pv.DataSet, field_name: str, expansion: float = 0.0
    ) -> tuple:
        """"""
        plot_config = self.plot_config

        cmap_choice = plot_config.cmap_mapping.get(field_name, "viridis")

        current_scalars = mesh.point_data[field_name] * plot_config.scale_mapping.get(
            field_name, 1.0
        )

        if field_name in ["viscosity", "strain_rate"]:
            current_scalars = np.log10(np.abs(np.maximum(current_scalars, 1e-30)))

        mesh[f"{field_name}_viz"] = current_scalars

        clim_actual = plot_config.clim_mapping.get(field_name, None)
        if clim_actual is None or clim_actual == "auto":
            finite_data = current_scalars[np.isfinite(current_scalars)]
            if finite_data.size > 0:
                clim_min = np.min(finite_data)
                clim_max = np.max(finite_data)
                if clim_min == clim_max:
                    clim_actual = (
                        clim_min - 0.1 * abs(clim_min) if clim_min != 0 else -0.1,
                        clim_max + 0.1 * abs(clim_max) if clim_max != 0 else 0.1,
                    )
                else:
                    if self._is_diverging_cmap(cmap_choice):
                        clim_max_abs = max(abs(clim_min), abs(clim_max))
                        clim_actual = (
                            (
                                -clim_max_abs - 0.1 * abs(clim_max_abs)
                                if clim_max_abs != 0
                                else -0.1
                            ),
                            (
                                clim_max_abs + 0.1 * abs(clim_max_abs)
                                if clim_max_abs != 0
                                else 0.1
                            ),
                        )
                    else:
                        clim_actual = (clim_min, clim_max)
            else:
                clim_actual = (0, 1)
            if self.verbosity >= 1:
                print(f" !! Warning: Using auto CLIM for '{field_name}': {clim_actual}")

        if expansion != 0.0:
            try:
                cmin, cmax = float(clim_actual[0]), float(clim_actual[1])
                crange = cmax - cmin
                clim_actual = (cmin - expansion * crange, cmax + expansion * crange)
            except (TypeError, ValueError) as e:
                if self.verbosity >= 1:
                    print(
                        f" !! Expansion skipped: CLIM values not numeric: {clim_actual} ({e})"
                    )

        modified_seq_target = ["stress_second_invariant", "X_field"]

        if self._is_diverging_cmap(cmap_choice):
            central_color = "#FFFFFF" if cmap_choice in ["RdGy"] else "#E5E5E5"
            final_cmap = self._modify_diverging_cmap(
                cmap_choice, central_color, plot_config.n_colors
            )
        elif (
            self._is_sequential_cmap(cmap_choice) and field_name in modified_seq_target
        ):
            final_cmap = self._modify_sequential_cmap(
                cmap_choice, "#E5E5E5", plot_config.n_colors
            )
        else:
            final_cmap = plt.get_cmap(cmap_choice, plot_config.n_colors)

        return final_cmap, clim_actual

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _compute_camera_settings(
        self, mesh, full_view: bool
    ) -> list[tuple[float, ...]]:
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
            bounds[2]
            + (bounds[3] - bounds[2]) * self.plot_config.scale_bar_position[1],
            0,
        ]
        end_point = [start_point[0] + scale_length, start_point[1], 0]

        line = pv.Line(start_point, end_point)

        plotter.add_mesh(
            line,
            color=self.plot_config.scale_bar_color,
            line_width=self.plot_config.scale_bar_thickness,
        )

        text_pos = [
            (start_point[0] + end_point[0]) / 2,
            start_point[1] + 0.01 * x_range,
            0,
        ]
        plotter.add_point_labels(
            [text_pos],
            [f"{scale_km:.0f} km"],
            font_size=self.plot_config.cbar_label_font_size,
            point_color=self.plot_config.scale_bar_color,
            text_color=self.plot_config.scale_bar_color,
            point_size=0,
            shape=None,
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def draw(self):
        """Main loop to process all models and generate visualizations."""
        warnings.filterwarnings("ignore", category=RuntimeWarning, module="pyvista")

        all_pvtu_files = self._get_pvtu_files()

        for model_id, result_files in all_pvtu_files.items():
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(f"==> Looking for results: {model_id}")

            if not result_files:
                if self.verbosity >= 1:
                    print(
                        f" !! Warning: no result files to process for model: {model_id}"
                    )
                continue

            if isinstance(self.plot_config.default_fig_dir, str):
                model_fig_dir = Path(self.plot_config.default_fig_dir)
            else:
                model_fig_dir = self.plot_config.default_fig_dir

            if not model_fig_dir.exists():
                model_fig_dir.mkdir(parents=True, exist_ok=True)

            timesteps_in_files = []
            for file in result_files:
                if isinstance(file, str):
                    file = Path(file)

                match = re.search(r"solution-(\d+)\.pvtu", file.name)
                if match:
                    timesteps_in_files.append(int(match.group(1)))
                else:
                    if self.verbosity >= 1:
                        print(
                            f" !! Warning: could not parse timestep from filename:\n"
                            f"    {file.name}"
                        )
                    timesteps_in_files.append(-1)

            for pvtu_path, tstep in zip(result_files, timesteps_in_files):
                if tstep == -1 or (self.tsteps and tstep not in self.tsteps):
                    continue

                if isinstance(pvtu_path, str):
                    pvtu_path = Path(pvtu_path)

                try:
                    mesh = pv.read(pvtu_path)
                    if mesh is None:
                        if self.verbosity >= 1:
                            print(
                                f" !! Warning: pyvista.read returned None for file:\n"
                                f"    {pvtu_path.name}"
                            )
                        continue
                except Exception as e:
                    if self.verbosity >= 1:
                        print(
                            f" !! Warning: failed to read mesh file {pvtu_path.name}:\n"
                            f" !! Error message: {e}"
                        )
                    continue

                print("    --------------------------------------------------")
                print(
                    f"    Drawing mesh plots @ timestep: {tstep} (File: {pvtu_path.name})"
                )
                print("    --------------------------------------------------")
                for field_key in self.plot_config.file_mapping.keys():
                    self.visualize_mesh(
                        mesh.copy(), field_key, tstep, model_id, model_fig_dir
                    )

                if self.plot_config.centerline_depth_plots:
                    print("    --------------------------------------------------")
                    print(
                        f"    Drawing centerline depth plots: @ timestep: {tstep} (File: {pvtu_path.name})"
                    )
                    print("    --------------------------------------------------")
                    for field_key in self.plot_config.file_mapping.keys():
                        self.visualize_centerline_depth(
                            mesh.copy(),
                            field_key,
                            tstep,
                            model_id,
                            model_fig_dir,
                            reaction_depth=132,
                        )

        pv.close_all()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def visualize_mesh(
        self,
        mesh: pv.DataSet,
        field_name: str,
        tstep: int,
        model_id: str,
        base_fig_dir: Path | str,
    ) -> None:
        """Visualizes a single field on the mesh using PyVista."""
        if isinstance(base_fig_dir, str):
            base_fig_dir = Path(base_fig_dir)

        plot_config = self.plot_config

        mesh = self._prepare_mesh_fields(mesh, field_name)
        mesh = self._rotate_mesh(mesh, tstep)

        if not self._check_mesh_field(mesh, field_name, model_id, tstep):
            return

        mapped_file_str = plot_config.file_mapping.get(
            field_name, field_name.replace("_", "-")
        )
        out_path = (
            base_fig_dir
            / f"{model_id.replace('_', '-')}-{mapped_file_str}-{str(tstep).zfill(4)}.png"
        )

        if out_path.exists():
            print(f" -- Found plot: {out_path.name}!")
            return

        print(f"--> {out_path.name}")

        time_myr = self._get_mesh_time_myr(mesh, model_id, tstep)

        plot_title = f"{int(time_myr):03} Ma"

        cmap, clim = self._configure_cmap(mesh, field_name)

        plotter: pv.Plotter = pv.Plotter(
            off_screen=True,
            window_size=plot_config.plotter_window_size,
            lighting=plot_config.plotter_lighting,
        )
        plotter.set_background(plot_config.plotter_background)  # type: ignore

        glyph_arrows = None
        if (
            field_name == "velocity" or plot_config.velocity_mapping[field_name]
        ) and "velocity" in mesh.point_data:
            if "velocity_mag" not in mesh.point_data:
                mesh["velocity_mag"] = np.linalg.norm(mesh["velocity"], axis=1)
            glyph_arrows = mesh.glyph(
                orient="velocity",  # type: ignore[arg-type]
                scale="velocity_mag",  # type: ignore[arg-type]
                factor=plot_config.velocity_glyph_factor,
                geom=pv.Arrow(),
                tolerance=0.05,
            )
        elif (
            field_name == "principal_stress_1"
            and "principal_stress_direction_1" in mesh.point_data
        ):
            glyph_arrows = mesh.glyph(
                orient="principal_stress_direction_1",  # type: ignore[arg-type]
                scale="principal_stress_1",  # type: ignore[arg-type]
                factor=plot_config.stress_glyph_factor,
                geom=pv.Line(),
            )
        elif (
            field_name == "principal_stress_2"
            and "principal_stress_direction_2" in mesh.point_data
        ):
            glyph_arrows = mesh.glyph(
                orient="principal_stress_direction_2",  # type: ignore[arg-type]
                scale="principal_stress_2",  # type: ignore[arg-type]
                factor=plot_config.stress_glyph_factor,
                geom=pv.Line(),
            )

        if glyph_arrows:
            plotter.add_mesh(
                glyph_arrows,
                color="black",
                line_width=plot_config.glyph_line_width,
                render_lines_as_tubes=False,
            )

        sargs = dict(
            title=plot_config.bar_mapping.get(field_name, field_name),
            vertical=plot_config.cbar_vertical,
            title_font_size=plot_config.cbar_title_font_size,
            label_font_size=plot_config.cbar_label_font_size,
            fmt=plot_config.fmt_mapping.get(field_name, "%.1f"),
            width=plot_config.cbar_width,
            n_labels=plot_config.cbar_n_labels,
            position_x=plot_config.cbar_position[0],
            position_y=plot_config.cbar_position[1],
        )

        plotter.add_mesh(
            mesh,
            scalars=f"{field_name}_viz",
            cmap=cmap,
            clim=clim,
            show_edges=plot_config.show_edges,
            edge_opacity=plot_config.edge_opacity,
            scalar_bar_args=sargs,
            nan_color="#FEFEFE",
        )

        camera_settings = self._compute_camera_settings(
            mesh, plot_config.camera_full_view
        )
        plotter.camera_position = camera_settings

        if plot_config.scale_bar_enabled:
            self._add_scale_bar(mesh, plotter)

        plotter.add_text(
            plot_title,
            font_size=plot_config.title_font_size,
            position=plot_config.title_position,  # type: ignore[arg-type]
            viewport=True,
        )

        if "depth" in mesh.point_data:
            for depth, width in zip(
                plot_config.depth_contour_depths_km,
                plot_config.depth_contour_line_widths,
            ):
                self._add_depth_contour(
                    plotter,
                    mesh,
                    depth,
                    line_width=width,
                    tolerance_km=plot_config.depth_contour_tolerance_km,
                )

        plotter.screenshot(out_path)
        plotter.close()

        try:
            img = Image.open(out_path)
            img.save(out_path, dpi=plot_config.screenshot_dpi)
        except Exception as e:
            print(f"PIL failed to resave {out_path} with new DPI: {e}")

        del mesh, plotter, glyph_arrows
        gc.collect()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def visualize_centerline_depth(
        self,
        mesh: pv.DataSet,
        field_name: str,
        tstep: int,
        model_id: str,
        base_fig_dir: Path | str,
        reaction_depth: float,
    ):
        """Plots a 1D vertical profile of a field through the center of the mesh in the x-direction."""
        if isinstance(base_fig_dir, str):
            base_fig_dir = Path(base_fig_dir)

        plot_config = self.plot_config

        mesh = self._prepare_mesh_fields(mesh, field_name)
        mesh = self._rotate_mesh(mesh, tstep)

        if not self._check_mesh_field(mesh, field_name, model_id, tstep):
            return

        mapped_file_str = plot_config.file_mapping.get(
            field_name, field_name.replace("_", "-")
        )
        out_path = (
            base_fig_dir
            / f"{model_id.replace('_', '-')}-centerline-{mapped_file_str}-{str(tstep).zfill(4)}.png"
        )

        if out_path.exists():
            print(f" -- Found plot: {out_path.name}!")
            return

        print(f"--> {out_path.name}")

        time_myr = self._get_mesh_time_myr(mesh, model_id, tstep)

        plot_title = f"{int(time_myr):03} Ma"

        cmap, clim = self._configure_cmap(mesh, field_name, 0.05)

        x_coords = mesh.points[:, 0]
        x_center = 0.5 * (x_coords.min() + x_coords.max())
        epsilon = (x_coords.max() - x_coords.min()) * 0.005
        center_mask = np.abs(x_coords - x_center) < epsilon

        if not np.any(center_mask):
            print(" !! No centerline points found!")
            return

        if (
            field_name == "velocity" or plot_config.velocity_mapping.get(field_name, False)
        ) and "velocity" in mesh.point_data:
            mesh["velocity"] = mesh["velocity"][:, 1]

        depths = mesh.point_data["depth"][center_mask] / 1e3
        values = mesh.point_data[field_name][
            center_mask
        ] * plot_config.scale_mapping.get(field_name, 1.0)

        sort_idx = np.argsort(depths)
        depths = depths[sort_idx]
        values = values[sort_idx]

        plt.rcParams.update(plot_config.centerline_plot_rcParams)
        fig, ax = plt.subplots(
            figsize=(
                plot_config.centerline_plot_fig_width,
                plot_config.centerline_plot_fig_height,
            )
        )

        ax.plot(values, depths, color=cmap(0.8))

        ax.set_xlabel(plot_config.bar_mapping.get(field_name, field_name))
        ax.set_ylabel("Depth (km)")
        ax.set_title(plot_title)
        ax.set_xlim(clim)
        ax.invert_yaxis()
        ax.grid(True, which="both", linewidth=0.5, color="#999999")
        ax.tick_params(axis="both", which="both", length=0)

        depth_max = depths.max()
        ax.axhspan(reaction_depth, depth_max, color="lightgrey", alpha=0.5, zorder=0)
        ax.axhline(reaction_depth, color="black", linestyle="--", linewidth=1, zorder=1)

        plt.tight_layout()
        plt.savefig(
            out_path, dpi=plot_config.centerline_plot_rcParams.get("figure.dpi", 300)
        )
        plt.close(fig)


#######################################################
## .3. Parse Arguments                           !!! ##
#######################################################
def parse_arguments() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Visualize simulation results.")
    parser.add_argument("--model-id", type=str, help="Simulation name")
    parser.add_argument("--timesteps", type=str, help="Timesteps to visualize")
    parser.add_argument("--in-dir", type=str, help="Simulation results dir")
    parser.add_argument("--out-fig-dir", type=str, help="Out fig dir")
    parser.add_argument("--visualization-type", type=str, help="Visualization type")

    return parser.parse_args()
