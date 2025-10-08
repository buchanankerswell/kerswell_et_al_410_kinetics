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
expand_range <- function(x, factor = 0.5) {
  r <- range(x, na.rm = TRUE)
  delta <- diff(r) * factor / 2
  c(r[1] - delta, r[2] + delta)
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
table_theme <- function() {
  theme(
    panel.ontop = TRUE,
    panel.grid.major = element_line(
      linewidth = 0.3,
      color = rgb(0.9, 0.9, 0.9, 0.8)
    ),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA),
    plot.margin = margin(50, 10, 5, 5),
    plot.title = element_text(hjust = 0.5),
    plot.tag.location = "panel",
    plot.tag = element_text(
      margin = margin(5, 0, 0, -5),
      hjust = 0,
      face = "bold",
      color = "white"
    ),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "inside",
    legend.position.inside = c(1.0, 1.0),
    legend.direction = "horizontal",
    legend.key.height = unit(1.0, "lines"),
    legend.key.width = unit(3.0, "lines"),
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
