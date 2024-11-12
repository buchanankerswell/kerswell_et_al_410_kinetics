#!/usr/bin/env Rscript

#######################################################
## .0. Load Libraries and Functions              !!! ##
#######################################################
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  if (length(match) > 0) {
    dirname(normalizePath(sub(needle, "", cmd_args[match])))
  } else if (!is.null(sys.frames()) && !is.null(sys.frame(1)$ofile)) {
    dirname(normalizePath(sys.frame(1)$ofile))
  } else {
    stop("Cannot determine script location.")
  }
}

this_dir <- get_script_dir()

source(file.path(this_dir, "utils.R"))

#######################################################
## .1. Visualize Adiabatic Profile               !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_adiabatic_profile <- function(profile_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df_profile <-
    read_burnman_profile(profile_path) |>
    mutate(
      pressure = pressure / 1e9,
      density = density / 1e3,
      thermal_expansivity = thermal_expansivity * 1e5,
      compressibility = compressibility * 1e12,
      specific_heat = specific_heat / 1e3
    ) |>
    select(-c(
      pressure,
      temperature,
      gravity,
      molar_gibbs,
      molar_entropy,
      molar_volume,
      seismic_vp,
      seismic_vs,
      seismic_dvp_dt,
      seismic_dvs_dt
    )) |>
    pivot_longer(
      cols = -depth,
      names_to = "variable",
      values_to = "value"
    )

  variable_order <- c(
    "pressure",
    "temperature",
    "density",
    "gravity",
    "thermal_expansivity",
    "compressibility",
    "specific_heat"
  )

  df_profile$variable <-
    factor(df_profile$variable, levels = variable_order)

  units <- c(
    "pressure" = "GPa",
    "temperature" = "K",
    "density" = "g/cm^3",
    "gravity" = "m/s^2",
    "thermal_expansivity" = "K%*%10^-5",
    "specific_heat" = "kJ/K~kg",
    "compressibility" = "Pa%*%10^-12"
  )

  labs <- c(
    "pressure" = "P",
    "temperature" = "T",
    "density" = "rho",
    "gravity" = "g",
    "thermal_expansivity" = "alpha",
    "specific_heat" = "Cp",
    "compressibility" = "beta"
  )

  custom_labeller <- function(variable) {
    paste0(labs[variable], " * ' [' * ", units[variable], " * ']'")
  }

  p <-
    df_profile |>
    ggplot(aes(x = value, y = depth / 1e6)) +
    geom_path(linewidth = 1.5) +
    facet_wrap(
      ~variable,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(variable = custom_labeller, .default = label_parsed)
    ) +
    scale_x_continuous(n.breaks = 5, expand = expansion(mult = 0.15)) +
    scale_y_reverse(expand = expansion(mult = 0.15)) +
    labs(x = NULL, y = bquote("Depth [" * km %*% 10^3 * "]")) +
    theme_bw(base_size = 36) +
    profile_theme() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold")
    )

  ggsave(
    out_path,
    plot = p,
    width = 20,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_driving_force_profile <- function(profile_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df_profile <-
    read_burnman_profile(profile_path) |>
    filter(depth <= 1000e3) |>
    mutate(
      depth = depth / 1e3,
      pressure = pressure / 1e9,
      molar_internal_energy_a = molar_internal_energy_a / 1e3,
      molar_internal_energy_b = molar_internal_energy_b / 1e3,
      delta_molar_internal_energy = delta_molar_internal_energy / 1e3,
      molar_gibbs_a = molar_gibbs_a / 1e3,
      molar_gibbs_b = molar_gibbs_b / 1e3,
      delta_molar_gibbs = delta_molar_gibbs / 1e3,
      molar_volume_a = molar_volume_a * 1e6,
      molar_volume_b = molar_volume_b * 1e6,
      delta_molar_volume = delta_molar_volume * 1e6
    ) |>
    select(
      depth,
      U_a = molar_internal_energy_a,
      U_b = molar_internal_energy_b,
      U_delta = delta_molar_internal_energy,
      G_a = molar_gibbs_a,
      G_b = molar_gibbs_b,
      G_delta = delta_molar_gibbs,
      S_a = molar_entropy_a,
      S_b = molar_entropy_b,
      S_delta = delta_molar_entropy,
      V_a = molar_volume_a,
      V_b = molar_volume_b,
      V_delta = delta_molar_volume
    ) |>
    pivot_longer(
      cols = -depth,
      names_to = c("property", "type"),
      names_sep = "_",
      values_to = "value"
    )

  units <- c(
    U = "kJ/mol",
    G = "kJ/mol",
    S = "J/K/mol",
    V = "cm^3/mol"
  )

  labs <- c(
    U = "U[eq]",
    G = "G[eq]",
    S = "S[eq]",
    V = "V[eq]"
  )

  labs_delta <- c(
    U = "Delta*U[eq]",
    G = "Delta*G[eq]",
    S = "Delta*S[eq]",
    V = "Delta*V[eq]"
  )

  custom_labeller <- function(variable) {
    paste0(labs[variable], " * ' [' * ", units[variable], " * ']'")
  }

  custom_labeller_delta <- function(variable) {
    paste0(labs_delta[variable], " * ' [' * ", units[variable], " * ']'")
  }

  df_phase <-
    filter(df_profile, type %in% c("a", "b")) |>
    mutate(type = ifelse(type == "a", "ol", "wad"))
  df_delta <-
    filter(df_profile, type == "delta") |>
    mutate(type = "(wad - ol)")

  p0 <-
    ggplot(df_phase, aes(x = value, y = depth, color = type)) +
    geom_path(linewidth = 1.5) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller, .default = label_parsed)
    ) +
    scale_y_reverse() +
    scale_color_brewer(palette = "Set1") +
    labs(x = NULL, y = "Depth [km]", color = NULL) +
    theme_bw(base_size = 36) +
    profile_theme() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold"),
      legend.position = "inside",
      legend.direction = "vertical",
      legend.position.inside = c(0.01, 0.15),
      legend.text = element_text(size = 28)
    )

  p1 <-
    ggplot(df_delta, aes(x = value, y = depth, color = type)) +
    geom_path(linewidth = 1.5) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller_delta, .default = label_parsed)
    ) +
    scale_y_reverse() +
    scale_color_manual(values = c("(wad - ol)" = "black")) +
    labs(x = NULL, y = bquote("Depth [km]"), color = NULL) +
    theme_bw(base_size = 36) +
    profile_theme() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold"),
      legend.position = "inside",
      legend.direction = "vertical",
      legend.position.inside = c(0.01, 0.90),
      legend.text = element_text(size = 28)
    )

  p <- p0 / p1

  ggsave(
    out_path,
    plot = p,
    width = 20,
    height = 13,
    dpi = 300,
    bg = "white"
  )
}
