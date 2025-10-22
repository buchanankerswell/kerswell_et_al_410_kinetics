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
to_snake_case <- function(s) {
  s |>
    tolower() |>
    str_replace_all(" ", "_")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_aspect_stats <- function(filepath) {
  file_lines <- readLines(filepath)
  header_lines <- file_lines[str_detect(file_lines, "^#")]
  col_names <-
    header_lines |>
    str_extract("(?<=: ).*") |>
    map_chr(to_snake_case)
  skip_lines <- length(header_lines)

  suppressWarnings({
    read_table(
      filepath,
      skip = skip_lines,
      col_names = col_names,
      show_col_types = FALSE
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_stokes_residuals <- function(filepath) {
  file_lines <- readLines(filepath)
  header_lines <- file_lines[str_detect(file_lines, "^#")]
  col_names <- c("time", "solveidx", "residual")
  skip_lines <- length(header_lines)

  suppressWarnings({
    df <- read_table(
      filepath,
      skip = skip_lines,
      col_names = col_names,
      show_col_types = FALSE
    )
  })

  seconds_in_a_year <- 365.25 * 24 * 60 * 60
  df$time <- df$time / seconds_in_a_year

  df |>
    group_by(time) |>
    summarize(
      min_residual = log10(min(residual)),
      max_residual = log10(max(residual)),
      mean_residual = log10(mean(residual))
    ) |>
    ungroup()
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_viscosity_profile <- function(filepath) {
  suppressWarnings({
    read_table(
      filepath,
      col_names = c("viscosity", "depth"),
      show_col_types = FALSE
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_lateral_profile <- function(filepath) {
  suppressWarnings({
    read_table(
      filepath,
      skip = 1,
      col_names = c("V", "depth"),
      show_col_types = FALSE
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_depth_average <- function(filepath) {
  file_lines <- readLines(filepath)
  header_lines <- file_lines[str_detect(file_lines, "^#")]
  col_names <- str_split_1(header_lines, "\\ +")
  col_names <- col_names[2:length(col_names)]
  skip_lines <- length(header_lines)

  suppressWarnings({
    df <- read_table(
      filepath,
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
read_displacements <- function(filepath) {
  suppressWarnings(df <- read_csv(filepath, show_col_types = FALSE))
  df |> mutate(displacement = displacement / 1e3, width = width / 1e3)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
profile_theme <- function() {
  theme(
    panel.grid.major = element_line(linewidth = 0.4, color = "white"),
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
theme_1 <- function() {
  theme(
    panel.grid.major = element_line(linewidth = 0.4, color = "white"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey90"),
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(),
    legend.justification = "left",
    legend.position = "inside",
    legend.position.inside = c(0.00, 0.85),
    legend.direction = "horizontal",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.box.margin = margin(2, 2, 2, 2),
    legend.margin = margin(),
    legend.title = element_text(vjust = 0, size = 14),
    legend.background = element_blank()
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_2 <- function() {
  theme(
    panel.grid.major = element_line(linewidth = 0.4, color = "white"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey90"),
    plot.margin = margin(5, 10, 5, 5),
    plot.title = element_text(hjust = 0.5),
    plot.tag.location = "panel",
    plot.tag.position = "topleft",
    plot.tag = element_text(size = 18, margin = margin(5, 0, 0, 0), hjust = 0, color = "black", face = "bold"),
    axis.ticks = element_blank(),
    legend.justification = "right",
    legend.position = "inside",
    legend.position.inside = c(0.93, 0.85),
    legend.direction = "horizontal",
    legend.key.height = unit(0.6, "lines"),
    legend.key.width = unit(1.2, "lines"),
    legend.ticks = element_line(color = "black", linewidth = 0.4),
    legend.ticks.length = unit(0.1, "lines"),
    legend.frame = element_rect(color = "black", linewidth = 0.4),
    legend.box.margin = margin(),
    legend.margin = margin(),
    legend.title = element_text(hjust = 1, vjust = 0, size = 12, margin = margin(0, 0, 2, 0)),
    legend.title.position = "top",
    legend.text = element_text(size = 11, margin = margin(2, 0, 0, 0)),
    legend.background = element_blank()
  )
}
