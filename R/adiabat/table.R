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
## .1. Visualize Material Table                  !!! ##
#######################################################
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
