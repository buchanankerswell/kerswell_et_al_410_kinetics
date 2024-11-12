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
## .1. Visualize Material Table                  !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_material_table <- function(profile_path, table_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df_profile <-
    read_burnman_profile(profile_path) |>
    mutate(
      pressure = pressure / 1e9,
      density = density / 1e3,
      thermal_expansivity = thermal_expansivity * 1e5,
      compressibility = compressibility * 1e12
    ) |>
    select(-c(seismic_vp, seismic_vs))

  df_tables <- read_material_table(table_path)

  props <- c(
    "density",
    "thermal_expansivity",
    "compressibility",
    "specific_heat"
  )

  units <- c(
    "density" = "g/cm^3",
    "thermal_expansivity" = "K%*%10^-5",
    "compressibility" = "Pa%*%10^-12",
    "specific_heat" = "J/K~kg"
  )

  labs <- c(
    "density" = "rho",
    "thermal_expansivity" = "alpha",
    "compressibility" = "beta",
    "specific_heat" = "Cp"
  )

  color_map <- c(
    "density" = "mako",
    "thermal_expansivity" = "magma",
    "compressibility" = "mako",
    "specific_heat" = "magma"
  )

  color_lims <- list(
    density = expand_range(df_profile$density),
    thermal_expansivity = expand_range(df_profile$thermal_expansivity),
    compressibility = expand_range(df_profile$compressibility),
    specific_heat = expand_range(df_profile$specific_heat)
  )

  color_lims <- list(
    density = expand_range_iqr(df_tables$density),
    thermal_expansivity = expand_range_iqr(df_tables$thermal_expansivity),
    compressibility = expand_range_iqr(df_tables$compressibility),
    specific_heat = expand_range_iqr(df_tables$specific_heat)
  )

  color_direction <- c(
    "density" = -1,
    "thermal_expansivity" = 1,
    "compressibility" = 1,
    "specific_heat" = 1
  )

  custom_labeller <- function(variable) {
    if (units[variable] != "") {
      label <- paste0(labs[variable], " * ' [' * ", units[variable], " * ']'")
    } else {
      label <- labs[variable]
    }
    parse(text = label)
  }

  suppressWarnings({
    p1 <-
      ggplot() +
      geom_raster(
        data = df_tables,
        aes(temperature, pressure, fill = get(props[1]))
      ) +
      geom_path(
        data = df_profile,
        aes(temperature, pressure),
        color = "white",
        linewidth = 1.4
      ) +
      scale_fill_viridis_c(
        option = color_map[props[1]],
        direction = color_direction[props[1]],
        limits = color_lims[[props[1]]],
        na.value = "grey90"
      ) +
      labs(
        x = "Temperature [K]",
        y = "Pressure [GPa]",
        fill = custom_labeller(props[1])
      ) +
      coord_cartesian(expand = FALSE) +
      theme_bw(base_size = 30) +
      table_theme() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      )

    p2 <-
      ggplot() +
      geom_raster(
        data = df_tables,
        aes(temperature, pressure, fill = get(props[2]))
      ) +
      geom_path(
        data = df_profile,
        aes(temperature, pressure),
        color = "white",
        linewidth = 1.4
      ) +
      scale_fill_viridis_c(
        option = color_map[props[2]],
        direction = color_direction[props[2]],
        limits = color_lims[[props[2]]],
        na.value = "grey90"
      ) +
      labs(
        x = "Temperature [K]",
        y = "Pressure [GPa]",
        fill = custom_labeller(props[2])
      ) +
      coord_cartesian(expand = FALSE) +
      theme_bw(base_size = 30) +
      table_theme() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      )

    p3 <-
      ggplot() +
      geom_raster(
        data = df_tables,
        aes(temperature, pressure, fill = get(props[3]))
      ) +
      geom_path(
        data = df_profile,
        aes(temperature, pressure),
        color = "white",
        linewidth = 1.4
      ) +
      scale_fill_viridis_c(
        option = color_map[props[3]],
        direction = color_direction[props[3]],
        limits = color_lims[[props[3]]],
        na.value = "grey90"
      ) +
      labs(
        x = "Temperature [K]",
        y = "Pressure [GPa]",
        fill = custom_labeller(props[3])
      ) +
      coord_cartesian(expand = FALSE) +
      theme_bw(base_size = 30) +
      table_theme()

    p4 <-
      ggplot() +
      geom_raster(
        data = df_tables,
        aes(temperature, pressure, fill = get(props[4]))
      ) +
      geom_path(
        data = df_profile,
        aes(temperature, pressure),
        color = "white",
        linewidth = 1.4
      ) +
      scale_fill_viridis_c(
        option = color_map[props[4]],
        direction = color_direction[props[4]],
        limits = color_lims[[props[4]]],
        na.value = "grey90"
      ) +
      labs(
        x = "Temperature [K]",
        y = "Pressure [GPa]",
        fill = custom_labeller(props[4])
      ) +
      coord_cartesian(expand = FALSE) +
      theme_bw(base_size = 30) +
      table_theme() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
      )

    p <- (p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "a")

    ggsave(
      out_path,
      plot = p,
      width = 13,
      height = 12,
      dpi = 300,
      bg = "white",
      create.dir = TRUE
    )
  })
}
