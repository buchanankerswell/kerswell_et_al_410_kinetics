#######################################################
## Visualize Adiabatic Reference Conditions      !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_burnman_profile <- function(filepath) {
  file_lines <- readLines(filepath)
  header_lines <- file_lines[str_detect(file_lines, "^#")]
  skip_lines <- length(header_lines)
  read_delim(
    filepath,
    delim = "\t",
    skip = skip_lines,
    show_col_types = FALSE
  ) |>
    mutate(across(everything(), ~ as.numeric(.)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_material_table <- function(filepath) {
  file_lines <- readLines(filepath)
  header_pattern <- "^T\\(K\\)\\s+P\\(bar\\)\\s+rho"
  header_line_indices <- str_which(file_lines, header_pattern)

  if (length(header_line_indices) == 1) {
    header_line_index <- header_line_indices[1]
    skip_lines <- header_line_index - 1
  } else if (length(header_line_indices) > 1) {
    stop("Multiple lines matched the header pattern!")
  } else {
    stop("Header line not found with the specified pattern!")
  }

  suppressWarnings({
    df <- read_table(filepath, skip = skip_lines, show_col_types = FALSE)
  })

  df |>
    rename(
      temperature = "T(K)",
      pressure = "P(bar)",
      density = "rho,kg/m3",
      thermal_expansivity = "alpha,1/K",
      compressibility = "beta,1/bar",
      specific_heat = "cp,J/K/kg",
      entropy = "s,J/K/kg"
    ) |>
    select(c(
      temperature,
      pressure,
      density,
      thermal_expansivity,
      compressibility,
      specific_heat,
      entropy
    )) |>
    mutate(across(everything(), ~ as.numeric(.))) |>
    mutate(
      pressure = pressure / 1e4,
      density = density / 1e3,
      thermal_expansivity = thermal_expansivity * 1e5,
      compressibility = compressibility * 1e7,
      specific_heat = specific_heat / 1e3,
      entropy = entropy / 1e3
    )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expand_range_iqr <- function(x, threshold = 5.0) {
  x <- x[is.finite(x)]
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - threshold * iqr
  upper <- q[2] + threshold * iqr
  inliers <- x[x >= lower & x <= upper]
  range(inliers, na.rm = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_profile <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey90"),
    panel.spacing.x = unit(1.5, "cm"),
    plot.margin = margin(5, 25, 5, 25),
    plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(),
    legend.justification = "left",
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.85),
    legend.direction = "horizontal",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.box.margin = margin(2, 2, 2, 2),
    legend.margin = margin(),
    legend.title = element_text(vjust = 0, size = 20),
    legend.background = element_blank()
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_table <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA),
    plot.margin = margin(50, 25, 5, 5),
    plot.title = element_text(hjust = 0.5),
    plot.tag.location = "panel",
    plot.tag = element_text(margin = margin(5, 0, 0, -5), hjust = 0, face = "bold", color = "white"),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "inside",
    legend.position.inside = c(1.0, 1.0),
    legend.direction = "horizontal",
    legend.key.height = unit(1.15, "lines"),
    legend.key.width = unit(2.0, "lines"),
    legend.ticks = element_line(color = "black", linewidth = 0.6),
    legend.ticks.length = unit(0.2, "lines"),
    legend.frame = element_rect(color = "black", linewidth = 0.6),
    legend.box.margin = margin(15, 5, 15, 5),
    legend.margin = margin(5, 5, 5, 5),
    legend.title = element_text(vjust = 1, size = 24),
    legend.text = element_text(size = 20, margin = margin(5, 0, 0, 0)),
    legend.background = element_blank()
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_thermodynamic_profile <- function(profile_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df_profile <- read_burnman_profile(profile_path) |>
    filter(pressure >= 10e9 & pressure <= 18e9) |>
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

  df_phase <- filter(df_profile, type %in% c("a", "b")) |>
    mutate(type = ifelse(type == "a", "ol", "wd"))
  df_delta <- filter(df_profile, type == "delta") |>
    mutate(type = "(wd - ol)")

  p0 <- ggplot(df_phase, aes(x = value, y = pressure, color = type)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller, .default = label_parsed)
    ) +
    scale_x_continuous(breaks = pretty_breaks(n = 4), guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.10)) +
    scale_y_reverse() +
    scale_color_brewer(palette = "Set1") +
    labs(x = NULL, y = "Pressure (GPa)", color = NULL) +
    theme_bw(base_size = 36) +
    theme_profile() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold"),
      legend.position = "inside",
      legend.direction = "vertical",
      legend.position.inside = c(0.01, 0.15),
      legend.text = element_text(size = 28)
    )

  p1 <- ggplot(df_delta, aes(x = value, y = pressure, color = type)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller_delta, .default = label_parsed)
    ) +
    scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.10)) +
    scale_y_reverse() +
    scale_color_manual(values = c("(wd - ol)" = "black")) +
    labs(x = NULL, y = bquote("Pressure (GPa)"), color = NULL) +
    theme_bw(base_size = 36) +
    theme_profile() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold"),
      legend.position = "inside",
      legend.direction = "vertical",
      legend.position.inside = c(0.01, 0.90),
      legend.text = element_text(size = 28)
    )

  p <- p0 / p1

  ggsave(out_path, plot = p, width = 15, height = 12, dpi = 300, bg = "white")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_material_profile <- function(profile_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df <- read_burnman_profile(profile_path) |>
    filter(pressure >= 10e9, pressure <= 18e9) |>
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
      type = recode(type, "a" = "ol", "b" = "wd")
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
      type = "(wd - ol)"
    )

  property_order <- c(
    "density",
    "gravity",
    "thermal_expansivity",
    "compressibility",
    "specific_heat"
  )

  df_phase$property <- factor(df_phase$property, levels = property_order)
  df_delta$property <- factor(df_delta$property, levels = property_order)

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

  p0 <- ggplot(df_phase, aes(x = value, y = pressure, color = type)) +
    geom_path(linewidth = 1.8) +
    facet_wrap(
      ~property,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(property = custom_labeller, .default = label_parsed)
    ) +
    scale_x_continuous(breaks = pretty_breaks(n = 4), guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.10)) +
    scale_y_reverse() +
    scale_color_brewer(palette = "Set1") +
    labs(x = NULL, y = "Pressure (GPa)", color = NULL) +
    theme_bw(base_size = 36) +
    theme_profile() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold"),
      legend.position = "inside",
      legend.direction = "vertical",
      legend.position.inside = c(0.01, 0.15),
      legend.text = element_text(size = 28)
    )

  ggsave(out_path, plot = p0, width = 20, height = 6.5, dpi = 300, bg = "white")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_material_table <- function(profile_path, table_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df_profile <- read_burnman_profile(profile_path) |>
    mutate(
      pressure = pressure / 1e9,
      density = density / 1e3,
      thermal_expansivity = thermal_expansivity * 1e5,
      compressibility = compressibility * 1e12
    ) |>
    select(-c(seismic_vp, seismic_vs)) |>
    filter(pressure >= 0 & pressure <= 25 & temperature >= 1553 & temperature <= 1923)

  df_tables <- read_material_table(table_path) |>
    filter(pressure >= 0 & pressure <= 25 & temperature >= 1553 & temperature <= 1923)

  props <- c("entropy", "density")

  units <- c(
    "density" = "g/cm^3",
    "thermal_expansivity" = "K%*%10^-5",
    "compressibility" = "Pa%*%10^-12",
    "specific_heat" = "KJ/K~kg",
    "entropy" = "kJ/K~kg"
  )

  labs <- c(
    "density" = "bar(rho)",
    "thermal_expansivity" = "bar(alpha)",
    "compressibility" = "bar(beta)",
    "specific_heat" = "bar(C[p])",
    "entropy" = "bar(S)"
  )

  color_map <- c(
    "density" = "mako",
    "thermal_expansivity" = "magma",
    "compressibility" = "mako",
    "specific_heat" = "magma",
    "entropy" = "mako"
  )

  color_lims <- list(
    density = expand_range_iqr(df_tables$density),
    thermal_expansivity = expand_range_iqr(df_tables$thermal_expansivity),
    compressibility = expand_range_iqr(df_tables$compressibility),
    specific_heat = expand_range_iqr(df_tables$specific_heat),
    entropy = expand_range_iqr(df_tables$entropy)
  )

  color_direction <- c(
    "density" = -1,
    "thermal_expansivity" = 1,
    "compressibility" = 1,
    "specific_heat" = 1,
    "entropy" = 1
  )

  custom_labeller <- function(variable) {
    if (units[variable] != "") {
      label <- paste0(labs[variable], " * ' (' * ", units[variable], " * ')'")
    } else {
      label <- labs[variable]
    }
    parse(text = label)
  }

  suppressWarnings({
    p1 <-
      ggplot() +
      geom_raster(data = df_tables, aes(temperature, pressure, fill = get(props[1]))) +
      geom_path(data = df_profile, aes(temperature, pressure), color = "white", linewidth = 1.4) +
      geom_path(data = data.frame(x = c(1705, 1809, 1809, 1705, 1705), y = c(10, 10, 18, 18, 10)), aes(x, y), color = "black", linewidth = 1.4) +
      scale_fill_viridis_c(
        option = color_map[props[1]],
        direction = color_direction[props[1]],
        limits = color_lims[[props[1]]],
        breaks = pretty_breaks(n = 3),
        na.value = "grey90"
      ) +
      labs(x = "Temperature (K)", y = "Pressure (GPa)", fill = custom_labeller(props[1])) +
      coord_cartesian(expand = FALSE) +
      theme_bw(base_size = 30) +
      theme_table()

    p2 <-
      ggplot() +
      geom_raster(data = df_tables, aes(temperature, pressure, fill = get(props[2]))) +
      geom_path(data = df_profile, aes(temperature, pressure), color = "white", linewidth = 1.4) +
      geom_path(data = data.frame(x = c(1705, 1809, 1809, 1705, 1705), y = c(10, 10, 18, 18, 10)), aes(x, y), color = "black", linewidth = 1.4) +
      scale_fill_viridis_c(
        option = color_map[props[2]],
        direction = color_direction[props[2]],
        limits = color_lims[[props[2]]],
        breaks = pretty_breaks(n = 4),
        na.value = "grey90"
      ) +
      labs(x = "Temperature (K)", y = "Pressure (GPa)", fill = custom_labeller(props[2])) +
      coord_cartesian(expand = FALSE) +
      theme_bw(base_size = 30) +
      theme_table() +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    p <- p1 + p2

    ggsave(out_path, plot = p, width = 10, height = 5.5, dpi = 300, bg = "white", create.dir = TRUE)
  })
}
