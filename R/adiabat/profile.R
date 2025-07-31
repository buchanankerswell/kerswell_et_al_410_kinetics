#!/usr/bin/env Rscript

#######################################################
## .0. Load Libraries and Functions              !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  if (length(match) > 0) {
    dirname(normalizePath(sub(needle, "", cmd_args[match])))
  } else if (!is.null(sys.frames()) && !is.null(sys.frame(1)$ofile)) {
    dirname(normalizePath(sys.frame(1)$ofile))
  } else {
    stop(" !! Error: cannot determine script location!")
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
    paste0(labs[variable], " * ' (' * ", units[variable], " * ')'")
  }

  p <-
    df_profile |>
    ggplot(aes(x = value, y = depth / 1e6)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~variable,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(variable = custom_labeller, .default = label_parsed)
    ) +
    scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.15)) +
    scale_y_reverse(expand = expansion(mult = 0.15)) +
    labs(x = NULL, y = bquote("Depth (" * km %*% 10^3 * ")")) +
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
    filter(pressure >= 10e9 & pressure <= 16e9) |>
    mutate(
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
      pressure,
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
      cols = -pressure,
      names_to = c("property", "type"),
      names_sep = "_",
      values_to = "value"
    )

  units <- c(
    G = "kJ/mol",
    S = "J/K/mol",
    V = "cm^3/mol"
  )

  labs <- c(
    G = "bar(G)",
    S = "bar(S)",
    V = "bar(V)"
  )

  labs_delta <- c(
    G = "Delta*bar(G)",
    S = "Delta*bar(S)",
    V = "Delta*bar(V)"
  )

  custom_labeller <- function(variable) {
    paste0(labs[variable], " * ' (' * ", units[variable], " * ')'")
  }

  custom_labeller_delta <- function(variable) {
    paste0(labs_delta[variable], " * ' (' * ", units[variable], " * ')'")
  }

  df_phase <-
    filter(df_profile, type %in% c("a", "b")) |>
    mutate(type = ifelse(type == "a", "ol", "wad"))
  df_delta <-
    filter(df_profile, type == "delta") |>
    mutate(type = "(wad - ol)")

  p0 <-
    ggplot(df_phase, aes(x = value, y = pressure, color = type)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller, .default = label_parsed)
    ) +
    scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.15)) +
    scale_y_reverse() +
    scale_color_brewer(palette = "Set1") +
    labs(x = NULL, y = "Pressure (GPa)", color = NULL) +
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
    ggplot(df_delta, aes(x = value, y = pressure, color = type)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller_delta, .default = label_parsed)
    ) +
    scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.15)) +
    scale_y_reverse() +
    scale_color_manual(values = c("(wad - ol)" = "black")) +
    labs(x = NULL, y = bquote("Pressure (GPa)"), color = NULL) +
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
    width = 15,
    height = 13,
    dpi = 300,
    bg = "white"
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_phase_transition_profile <- function(profile_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df <-
    read_burnman_profile(profile_path) |>
    filter(pressure >= 10e9, pressure <= 16e9) |>
    select(-c(
      molar_internal_energy_a,
      molar_internal_energy_b,
      delta_molar_internal_energy,
      molar_gibbs_a,
      molar_gibbs_b,
      delta_molar_gibbs,
      molar_entropy_a,
      molar_entropy_b,
      delta_molar_entropy,
      molar_volume_a,
      molar_volume_b,
      delta_molar_volume,
      pressure_wave_velocity_a,
      pressure_wave_velocity_b,
      pressure_wave_velocity_T_derivative_a,
      pressure_wave_velocity_T_derivative_b,
      shear_wave_velocity_a,
      shear_wave_velocity_b,
      shear_wave_velocity_T_derivative_a,
      shear_wave_velocity_T_derivative_b
    ))

  df_phase <- df |>
    select(pressure, matches("_[ab]$")) |>
    pivot_longer(
      cols = -pressure,
      names_to = c("property", "type"),
      names_sep = "_(?=[ab]$)",
      values_to = "value"
    ) |>
    mutate(
      pressure = pressure / 1e9,
      value = case_when(
        grepl("density", property) ~ value / 1e3,
        grepl("thermal_expansivity", property) ~ value * 1e5,
        grepl("specific_heat", property) ~ value / 1e3,
        grepl("compressibility", property) ~ value * 1e12,
        TRUE ~ value
      ),
      type = recode(type, "a" = "ol", "b" = "wad")
    )

  df_delta <- df |>
    select(pressure, starts_with("delta_")) |>
    pivot_longer(
      cols = -pressure,
      names_to = "property",
      values_to = "value"
    ) |>
    mutate(
      pressure = pressure / 1e9,
      property = sub("^delta_", "", property),
      value = case_when(
        property == "density" ~ value / 1e3,
        property == "thermal_expansivity" ~ value * 1e5,
        property == "specific_heat" ~ value / 1e3,
        property == "compressibility" ~ value * 1e12,
        TRUE ~ value
      ),
      type = "(wad - ol)"
    )

  property_order <- c(
    "density",
    "gravity",
    "thermal_expansivity",
    "compressibility",
    "specific_heat"
  )

  df_phase$property <-
    factor(df_phase$property, levels = property_order)
  df_delta$property <-
    factor(df_delta$property, levels = property_order)

  units <- c(
    "density" = "g/cm^3",
    "thermal_expansivity" = "K%*%10^-5",
    "specific_heat" = "kJ/K~kg",
    "compressibility" = "Pa%*%10^-12"
  )

  labs <- c(
    "density" = "bar(rho)",
    "thermal_expansivity" = "bar(alpha)",
    "specific_heat" = "bar(C)[p]",
    "compressibility" = "bar(beta)"
  )

  labs_delta <- c(
    "density" = "Delta*bar(rho)",
    "thermal_expansivity" = "Delta*bar(alpha)",
    "specific_heat" = "Delta*bar(C)[p]",
    "compressibility" = "Delta*bar(beta)"
  )

  custom_labeller <- function(variable) {
    paste0(labs[variable], " * ' (' * ", units[variable], " * ')'")
  }

  custom_labeller_delta <- function(variable) {
    paste0(labs_delta[variable], " * ' (' * ", units[variable], " * ')'")
  }

  p0 <-
    ggplot(df_phase, aes(x = value, y = pressure, color = type)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller, .default = label_parsed)
    ) +
    scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.15)) +
    scale_y_reverse() +
    scale_color_brewer(palette = "Set1") +
    labs(x = NULL, y = "Pressure (GPa)", color = NULL) +
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
    ggplot(df_delta, aes(x = value, y = pressure, color = type)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller_delta, .default = label_parsed)
    ) +
    scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.15)) +
    scale_y_reverse() +
    scale_color_manual(values = c("(wad - ol)" = "black")) +
    labs(x = NULL, y = bquote("Pressure (GPa)"), color = NULL) +
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

  p <- p0
  # p <- p0 / p1

  ggsave(
    out_path,
    plot = p,
    width = 20,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}
