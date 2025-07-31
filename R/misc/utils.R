#!/usr/bin/env Rscript

#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(zoo)
  library(tools)
})

#######################################################
## .1. Utility Functions                         !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ensure_output_dir <- function(out_path) {
  parent_dir <- dirname(out_path)
  if (!dir.exists(parent_dir)) {
    dir.create(parent_dir)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_exists <- function(out_path) {
  if (file.exists(out_path)) {
    cat(" -- Found plot: ", basename(out_path), "!\n", sep = "")
    return(TRUE)
  }

  ensure_output_dir(out_path)

  cat(" -> ", basename(out_path), "\n", sep = "")
  FALSE
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calc_viscosity <- function(pressure, temperature, eps, n, m, d, ea, va, a) {
  gas_const <- 8.3144
  shear_modulus <- 8e10
  burgers_vector <- 5e-10

  0.5 * shear_modulus * a^(-1 / n) *
    (d / burgers_vector)^(m / n) * eps^((1 - n) / n) *
    exp((ea + (va * pressure)) / (n * gas_const * temperature))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calc_stress <- function(pressure, temperature, eps, n, m, d, ea, va, a) {
  gas_const <- 8.3144
  shear_modulus <- 8e10
  burgers_vector <- 5e-10

  shear_modulus * a^(-1 / n) * (d / burgers_vector)^(m / n) * eps^(1 / n) *
    exp((ea + (va * pressure)) / (n * gas_const * temperature))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_depth_average <- function(file_path) {
  file_lines <- readLines(file_path)
  header_lines <- file_lines[str_detect(file_lines, "^#")]
  col_names <- str_split_1(header_lines, "\\ +")
  col_names <- col_names[2:length(col_names)]
  skip_lines <- length(header_lines)

  suppressWarnings({
    df <- read_table(
      file_path,
      skip = skip_lines,
      col_names = col_names,
      show_col_types = FALSE
    )
  })

  seconds_in_a_year <- 365.25 * 24 * 60 * 60

  df |> mutate(
    time = time / seconds_in_a_year,
    sinking_velocity = sinking_velocity * 1e2,
    vertical_heat_flux = vertical_heat_flux * 1e3,
    vertical_mass_flux = vertical_mass_flux * seconds_in_a_year * 1e-1
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
profile_theme <- function() {
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "white"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey90"),
    plot.margin = margin(15, 15, 15, 15),
    plot.title = element_text(hjust = 0.5),
    plot.tag.location = "panel",
    plot.tag = element_text(margin = margin(2, 2, 2, -4), hjust = 0),
    axis.ticks = element_blank(),
    legend.justification = "left",
    legend.position = "bottom",
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
deformation_theme <- function() {
  theme(
    panel.grid.major = element_line(linewidth = 0.3, color = "white"),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "grey90"),
    plot.margin = margin(15, 15, 15, 15),
    plot.title = element_text(hjust = 0.5),
    plot.tag.location = "panel",
    plot.tag = element_text(margin = margin(5, 2, 2, -5), hjust = 0),
    axis.ticks = element_blank(),
    legend.justification = "left",
    legend.position = "inside",
    legend.position.inside = c(0.70, 0.64),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.box.margin = margin(2, 2, 2, 2),
    legend.margin = margin(),
    legend.title = element_text(
      vjust = 0,
      size = 28,
      margin = margin(0, 0, 5, 0)
    ),
    legend.background = element_rect(fill = "grey90")
  )
}
