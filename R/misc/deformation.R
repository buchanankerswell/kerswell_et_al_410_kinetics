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
## .1. Visualize Deformation Map                 !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_deformation_map <- function(out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  temperature <- seq(273, 3773)
  strain_rate <- c(1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12)

  df_grid <-
    expand_grid(temperature, strain_rate) |>
    mutate(pressure = 0)

  d_diff_dry <-
    df_grid |>
    rowwise() |>
    mutate(
      viscosity = calc_viscosity(
        pressure,
        temperature,
        strain_rate,
        n = 1,
        m = 2.5,
        d = 1e-3,
        ea = 300e3,
        va = 6e-6,
        a = 8.7e15
      ),
      stress = calc_stress(
        pressure,
        temperature,
        strain_rate,
        n = 1,
        m = 2.5,
        d = 1e-3,
        ea = 300e3,
        va = 6e-6,
        a = 8.7e15
      )
    ) |>
    ungroup()

  d_diff_wet <-
    df_grid |>
    rowwise() |>
    mutate(
      viscosity = calc_viscosity(
        pressure,
        temperature,
        strain_rate,
        n = 1,
        m = 2.5,
        d = 1e-3,
        ea = 240e3,
        va = 5e-6, a = 5.3e15
      ),
      stress = calc_stress(
        pressure,
        temperature,
        strain_rate,
        n = 1,
        m = 2.5,
        d = 1e-3,
        ea = 240e3,
        va = 5e-6, a = 5.3e15
      )
    ) |>
    ungroup()

  d_disl_dry <-
    df_grid |>
    rowwise() |>
    mutate(
      viscosity = calc_viscosity(
        pressure,
        temperature,
        strain_rate,
        n = 3.5,
        m = 0,
        d = 1e-3,
        ea = 540e3,
        va = 15e-6,
        a = 3.5e22
      ),
      stress = calc_stress(
        pressure,
        temperature,
        strain_rate,
        n = 3.5,
        m = 0,
        d = 1e-3,
        ea = 540e3,
        va = 15e-6,
        a = 3.5e22
      )
    ) |>
    ungroup()

  d_disl_wet <-
    df_grid |>
    rowwise() |>
    mutate(
      viscosity = calc_viscosity(
        pressure,
        temperature,
        strain_rate,
        n = 3.5,
        m = 0,
        d = 1e-3,
        ea = 430e3,
        va = 10e-6,
        a = 2.0e18
      ),
      stress = calc_stress(
        pressure,
        temperature,
        strain_rate,
        n = 3.5,
        m = 0,
        d = 1e-3,
        ea = 430e3,
        va = 10e-6,
        a = 2.0e18
      )
    ) |>
    ungroup()

  d_comb_dry <-
    tibble(
      temperature = d_diff_dry$temperature,
      strain_rate = d_diff_dry$strain_rate,
      viscosity_diff = d_diff_dry$viscosity,
      viscosity_disl = d_disl_dry$viscosity,
      stress_diff = d_diff_dry$stress,
      stress_disl = d_disl_dry$stress
    )

  d_comb_wet <-
    tibble(
      temperature = d_diff_wet$temperature,
      strain_rate = d_diff_wet$strain_rate,
      viscosity_diff = d_diff_wet$viscosity,
      viscosity_disl = d_disl_wet$viscosity,
      stress_diff = d_diff_wet$stress,
      stress_disl = d_disl_wet$stress
    )

  d_eff_dry <-
    d_comb_dry |>
    rowwise() |>
    mutate(
      viscosity = min(viscosity_diff, viscosity_disl),
      stress = min(stress_diff, stress_disl),
      min_viscosity_type = if_else(
        viscosity == viscosity_diff, "diff",
        "disl"
      ),
      min_stress_type = if_else(
        stress == stress_diff, "diff",
        "disl"
      )
    ) |>
    ungroup()

  d_eff_wet <-
    d_comb_wet |>
    rowwise() |>
    mutate(
      viscosity = min(viscosity_diff, viscosity_disl),
      stress = min(stress_diff, stress_disl),
      min_viscosity_type = if_else(
        viscosity == viscosity_diff, "diff",
        "disl"
      ),
      min_stress_type = if_else(
        stress == stress_diff, "diff",
        "disl"
      )
    ) |>
    ungroup()

  suppressWarnings({
    p1 <-
      ggplot(mapping = aes(
        x = temperature,
        y = log10(stress),
        group = factor(strain_rate),
        linewidth = factor(strain_rate)
      )) +
      geom_path(data = d_eff_dry, linetype = 1) +
      geom_path(
        data = filter(d_eff_dry, strain_rate == 1e-15),
        linetype = 1,
        color = "deeppink"
      ) +
      geom_point(
        data = filter(d_eff_dry, strain_rate == 1e-15, temperature == 1600),
        size = 5,
        shape = 18
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = "Dry",
        hjust = 1.2,
        vjust = 1.4,
        size = 8
      ) +
      annotate(
        "text",
        x = 1800,
        y = 7,
        label = "Dislocation\nCreep",
        size = 8
      ) +
      annotate(
        "text",
        x = 1100,
        y = 3.5,
        label = "Diffusion\nCreep",
        size = 8
      ) +
      scale_x_continuous(
        limits = c(1000, 2000),
        breaks = seq(1000, 2000, 200),
        expand = expansion(mult = 0.15)
      ) +
      scale_y_continuous(
        limits = c(3, 8),
        expand = expansion(mult = 0.15)
      ) +
      scale_linewidth_discrete(range = c(0.3, 1.3)) +
      labs(
        x = "Temperature [K]",
        y = bquote("Log" ~ sigma ~ "[Pa]"),
        linewidth = bquote(dot(epsilon) ~ "[" * s^-1 * "]")
      ) +
      theme_bw(base_size = 30) +
      deformation_theme() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      )

    p2 <-
      ggplot(mapping = aes(
        x = temperature,
        y = log10(viscosity),
        group = factor(strain_rate),
        linewidth = factor(strain_rate)
      )) +
      geom_path(data = d_eff_dry, linetype = 1) +
      geom_path(
        data = filter(d_eff_dry, strain_rate == 1e-15),
        linetype = 1,
        color = "deeppink"
      ) +
      geom_point(
        data = filter(d_eff_dry, strain_rate == 1e-15, temperature == 1600),
        size = 5,
        shape = 18
      ) +
      annotate(
        "text",
        x = -Inf,
        y = -Inf,
        label = "Dry",
        hjust = -0.2,
        vjust = -0.4,
        size = 8
      ) +
      scale_x_continuous(
        limits = c(1000, 2000),
        breaks = seq(1000, 2000, 200),
        expand = expansion(mult = 0.15)
      ) +
      scale_y_continuous(
        limits = c(16, 25),
        expand = expansion(mult = 0.15),
        labels = scales::label_number(accuracy = 1)
      ) +
      scale_linewidth_discrete(range = c(0.3, 1.3)) +
      labs(
        x = "Temperature [K]",
        y = bquote("Log" ~ sigma ~ "[Pa]"),
        linewidth = bquote(dot(epsilon) ~ "[" * s^-1 * "]")
      ) +
      theme_bw(base_size = 30) +
      deformation_theme() +
      theme(legend.position = "none")

    p3 <-
      ggplot(mapping = aes(
        x = temperature,
        y = log10(stress),
        group = factor(strain_rate),
        linewidth = factor(strain_rate)
      )) +
      geom_path(data = d_eff_wet, linetype = 1) +
      geom_path(
        data = filter(d_eff_wet, strain_rate == 1e-15),
        linetype = 1,
        color = "deeppink"
      ) +
      geom_point(
        data = filter(d_eff_wet, strain_rate == 1e-15, temperature == 1600),
        size = 5,
        shape = 18
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = "Wet",
        hjust = 1.2,
        vjust = 1.4,
        size = 8
      ) +
      scale_x_continuous(
        limits = c(1000, 2000),
        breaks = seq(1000, 2000, 200),
        expand = expansion(mult = 0.15)
      ) +
      scale_y_continuous(
        limits = c(3, 8),
        expand = expansion(mult = 0.15)
      ) +
      scale_linewidth_discrete(range = c(0.3, 1.3)) +
      labs(
        x = "Temperature [K]",
        y = bquote("Log" ~ sigma ~ "[Pa]"),
        linewidth = bquote(dot(epsilon) ~ "[" * s^-1 * "]")
      ) +
      theme_bw(base_size = 30) +
      deformation_theme() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      )

    p4 <-
      ggplot(mapping = aes(
        x = temperature,
        y = log10(viscosity),
        group = factor(strain_rate),
        linewidth = factor(strain_rate)
      )) +
      geom_path(data = d_eff_wet, linetype = 1) +
      geom_path(
        data = filter(d_eff_wet, strain_rate == 1e-15),
        linetype = 1,
        color = "deeppink"
      ) +
      geom_point(
        data = filter(d_eff_wet, strain_rate == 1e-15, temperature == 1600),
        size = 5,
        shape = 18
      ) +
      annotate(
        "text",
        x = -Inf,
        y = -Inf,
        label = "Wet",
        hjust = -0.2,
        vjust = -0.4,
        size = 8
      ) +
      scale_x_continuous(
        limits = c(1000, 2000),
        breaks = seq(1000, 2000, 200),
        expand = expansion(mult = 0.15)
      ) +
      scale_y_continuous(
        limits = c(16, 25),
        expand = expansion(mult = 0.15),
        labels = scales::label_number(accuracy = 1)
      ) +
      scale_linewidth_discrete(range = c(0.3, 1.3)) +
      labs(
        x = "Temperature [K]",
        y = bquote("Log" ~ sigma ~ "[Pa]"),
        linewidth = bquote(dot(epsilon) ~ "[" * s^-1 * "]")
      ) +
      theme_bw(base_size = 30) +
      deformation_theme() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      )

    p <- (p1 + p3) / (p2 + p4) + plot_annotation(tag_levels = "a")

    ggsave(
      out_path,
      plot = p,
      width = 13,
      height = 11,
      dpi = 300,
      bg = "white"
    )
  })
}
