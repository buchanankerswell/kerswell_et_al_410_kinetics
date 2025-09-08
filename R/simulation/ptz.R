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
## .1. Visualize displacement                      !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_displacement <- function(in_path, out_path1, out_path2, out_path3) {
  if (plot_exists(out_path1) && plot_exists(out_path2)) {
    return(invisible())
  }

  df <-
    read_displacements(in_path) |>
    filter(timestep == 100) |>
    mutate(model_id = str_replace_all(model_id, "_", "-"))
  df_slab <- df |>
    filter(str_detect(model_id, "slab")) %>%
    mutate(model_id = str_replace_all(model_id, "slab-", "")) |>
    mutate(model_id = str_replace_all(model_id, "-", " "))
  df_plume <- df |>
    filter(str_detect(model_id, "plume")) %>%
    mutate(model_id = str_replace_all(model_id, "plume-", "")) |>
    mutate(model_id = str_replace_all(model_id, "-", " ")) |>
    mutate(width = abs(width))

  df_slab_no_anomalies <- df_slab |> filter(displacement >= -85)

  linear_model_plume <- lm(width ~ displacement, data = df_plume)
  log_model_plume <- lm(log10(width) ~ log10(max_reaction_rate), data = df_plume)
  y_shift_plume <- abs(min(df_plume$displacement)) + 1
  df_plume_shifted <- df_plume |> mutate(displacement_shifted = displacement + y_shift_plume)
  log_model_plume2 <- lm(log10(displacement + y_shift_plume) ~ log10(max_reaction_rate), data = df_plume_shifted)
  quadratic_model_slab <- lm(width ~ displacement + I(displacement^2), data = df_slab)
  quadratic_model_slab_no_anomalies <- lm(width ~ displacement + I(displacement^2), data = df_slab_no_anomalies)
  y_shift_slab <- abs(min(df_slab$displacement)) + 1
  df_slab_shifted <- df_slab |> mutate(displacement_shifted = displacement + y_shift_slab)
  df_slab_shifted_no_anomalies <- df_slab_no_anomalies |> mutate(displacement_shifted = displacement + y_shift_slab)
  log_model_slab <- lm(log10(displacement + y_shift_slab) ~ log10(max_reaction_rate), data = df_slab_shifted)
  log_model_slab_no_anomalies <- lm(log10(displacement + y_shift_slab) ~ log10(max_reaction_rate), data = df_slab_shifted_no_anomalies)
  log_quadratic_model_slab <- lm(log10(width) ~ poly(log10(max_reaction_rate), 2, raw = TRUE), data = df_slab)
  log_quadratic_model_slab_no_anomalies <- lm(log10(width) ~ poly(log10(max_reaction_rate), 2, raw = TRUE), data = df_slab_no_anomalies)

  width_range_plume <- seq(min(df_plume$width), max(df_plume$width), length.out = 100)
  displacement_range_plume <- seq(min(df_plume$displacement), max(df_plume$displacement), length.out = 100)
  max_reaction_rate_range_plume <- seq(min(df_plume$max_reaction_rate), max(df_plume$max_reaction_rate), length.out = 100)

  width_range_slab <- seq(-90, max(df_slab$width), length.out = 100)
  displacement_range_slab <- seq(-90, max(df_slab$displacement), length.out = 100)
  max_reaction_rate_range_slab <- seq(min(df_slab$max_reaction_rate), max(df_slab$max_reaction_rate), length.out = 100)

  linear_plume_predictions <- data.frame(
    displacement = displacement_range_plume,
    width = predict(linear_model_plume, newdata = data.frame(displacement = displacement_range_plume))
  )
  log_plume_predictions <- data.frame(
    max_reaction_rate = max_reaction_rate_range_plume,
    width = 10^predict(log_model_plume, newdata = data.frame(max_reaction_rate = max_reaction_rate_range_plume))
  )
  slab_predictions <- data.frame(
    displacement = displacement_range_slab,
    width = predict(quadratic_model_slab, newdata = data.frame(displacement = displacement_range_slab))
  )
  slab_predictions_no_anomalies <- data.frame(
    displacement = displacement_range_slab,
    width = predict(quadratic_model_slab_no_anomalies, newdata = data.frame(displacement = displacement_range_slab))
  )
  log_quadratic_slab_predictions <- data.frame(
    max_reaction_rate = max_reaction_rate_range_slab,
    width = 10^predict(log_quadratic_model_slab, newdata = data.frame(max_reaction_rate = max_reaction_rate_range_slab))
  )
  log_quadratic_slab_predictions_no_anomalies <- data.frame(
    max_reaction_rate = max_reaction_rate_range_slab,
    width = 10^predict(log_quadratic_model_slab_no_anomalies, newdata = data.frame(max_reaction_rate = max_reaction_rate_range_slab))
  )

  linear_plume_coef <- coef(linear_model_plume)
  linear_plume_r2 <- summary(linear_model_plume)$r.squared

  log_plume_coef <- coef(log_model_plume)
  log_plume_r2 <- summary(log_model_plume)$r.squared

  log_plume_coef2 <- coef(log_model_plume2)
  log_plume_r22 <- summary(log_model_plume2)$r.squared

  quadratic_slab_coef <- coef(quadratic_model_slab_no_anomalies)
  quadratic_slab_r2 <- summary(quadratic_model_slab_no_anomalies)$r.squared

  log_slab_coef <- coef(log_model_slab)
  log_slab_r2 <- summary(log_model_slab)$r.squared

  log_slab_coef_no_anomalies <- coef(log_model_slab_no_anomalies)
  log_slab_r2_no_anomalies <- summary(log_model_slab_no_anomalies)$r.squared

  log_quadratic_slab_coef <- coef(log_quadratic_model_slab_no_anomalies)
  log_quadratic_slab_r2 <- summary(log_quadratic_model_slab_no_anomalies)$r.squared

  log_plume_predictions2 <- data.frame(
    max_reaction_rate = max_reaction_rate_range_plume,
    displacement = 10^log_plume_coef2[1] * max_reaction_rate_range_plume^log_plume_coef2[2] - y_shift_plume
  )
  log_slab_predictions <- data.frame(
    max_reaction_rate = max_reaction_rate_range_slab,
    displacement = 10^log_slab_coef[1] * max_reaction_rate_range_slab^log_slab_coef[2] - y_shift_slab
  )
  log_slab_predictions_no_anomalies <- data.frame(
    max_reaction_rate = max_reaction_rate_range_slab,
    displacement = 10^log_slab_coef_no_anomalies[1] * max_reaction_rate_range_slab^log_slab_coef_no_anomalies[2] - y_shift_slab
  )


  linear_plume_equation <- paste0(
    "R² = ", round(linear_plume_r2, 3),
    "\ny = ", round(linear_plume_coef[1], 1), " + ", round(linear_plume_coef[2], 2), "x"
  )
  log_plume_equation <- paste0(
    "R² = ", round(log_plume_r2, 3),
    "\ny = ", round(10^log_plume_coef[1], 1), " * x^", round(log_plume_coef[2], 2)
  )
  log_plume_equation2 <- paste0(
    "R² = ", round(log_plume_r22, 3),
    "\ny = ", round(10^log_plume_coef2[1], 1), " * x^", round(log_plume_coef2[2], 2),
    " + ", round(-y_shift_plume, 1)
  )
  log_slab_equation <- paste0(
    "R² = ", round(log_slab_r2_no_anomalies, 3),
    "\ny = ", round(10^log_slab_coef_no_anomalies[1], 1), " * x^", round(log_slab_coef_no_anomalies[2], 2),
    " + ", round(-y_shift_slab, 1)
  )
  quadratic_slab_equation <- paste0(
    "R² = ", round(quadratic_slab_r2, 3),
    "\ny = ", round(quadratic_slab_coef[1], 1), " + ", round(quadratic_slab_coef[2], 2), "x + ", round(quadratic_slab_coef[3], 2), "x²"
  )
  log_quadratic_slab_equation <- paste0(
    "R² = ", round(log_quadratic_slab_r2, 3),
    "\ny = ", round(10^log_quadratic_slab_coef[1], 2), " * x^", round(log_quadratic_slab_coef[2], 2), " *\n",
    "10^(", round(log_quadratic_slab_coef[3], 3), " * log(x)²)"
  )

  p0 <-
    df_plume |>
    ggplot(aes(x = displacement, y = width)) +
    geom_line(
      data = linear_plume_predictions, aes(x = displacement, y = width),
      color = "black", linewidth = 0.8, inherit.aes = FALSE
    ) +
    geom_point(color = "black") +
    annotate("text",
      x = 74, y = 5,
      label = linear_plume_equation,
      hjust = 1, vjust = 0,
      size = 4
    ) +
    labs(
      x = "PTZ Displacement (km)",
      y = "PTZ Width (km)"
    ) +
    theme_bw(base_size = 14) +
    theme_2()


  p1 <-
    filter(df_slab, displacement > -85) |>
    ggplot(aes(x = displacement, y = width)) +
    geom_line(
      data = slab_predictions, aes(x = displacement, y = width),
      color = "black", linetype = 2, linewidth = 0.5, inherit.aes = FALSE
    ) +
    geom_line(
      data = slab_predictions_no_anomalies, aes(x = displacement, y = width),
      color = "black", linewidth = 0.8, inherit.aes = FALSE
    ) +
    geom_point(data = filter(df_slab, displacement <= -85), color = "red") +
    geom_point(color = "black") +
    annotate("text",
      x = -112, y = 4.9,
      label = quadratic_slab_equation,
      hjust = 0, vjust = 0,
      size = 4
    ) +
    labs(
      x = "PTZ Displacement (km)",
      y = "PTZ Width (km)"
    ) +
    theme_bw(base_size = 14) +
    theme_2()

  p2 <-
    df_plume |>
    ggplot(aes(x = max_reaction_rate, y = width)) +
    geom_line(
      data = log_plume_predictions, aes(x = max_reaction_rate, y = width),
      color = "black", linewidth = 0.8, inherit.aes = FALSE
    ) +
    geom_point(color = "black") +
    annotate("text",
      x = ((max(df_plume$max_reaction_rate) - min(df_plume$max_reaction_rate)) / 2) + 1.5,
      y = (max(df_plume$width) - min(df_plume$width)) / 2,
      label = log_plume_equation,
      hjust = 0.5, vjust = 0.5,
      size = 4
    ) +
    labs(
      x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      y = "PTZ Width (km)"
    ) +
    theme_bw(base_size = 14) +
    theme_2() +
    theme(
      plot.tag.position = "topright",
      plot.tag = element_text(margin = margin(5, 10, 0, 0))
    )

  p3 <-
    filter(df_slab, displacement > -85) |>
    ggplot(aes(x = max_reaction_rate, y = width)) +
    geom_line(
      data = log_quadratic_slab_predictions, aes(x = max_reaction_rate, y = width),
      color = "black", linetype = 2, linewidth = 0.5, inherit.aes = FALSE
    ) +
    geom_line(
      data = log_quadratic_slab_predictions_no_anomalies, aes(x = max_reaction_rate, y = width),
      color = "black", linewidth = 0.8, inherit.aes = FALSE
    ) +
    geom_point(data = filter(df_slab, displacement <= -85), color = "red") +
    geom_point(color = "black") +
    annotate("text",
      x = ((max(df_slab$max_reaction_rate) - min(df_slab$max_reaction_rate)) / 2) + 0.2,
      y = (max(df_plume$width) - min(df_plume$width)) / 2,
      label = log_quadratic_slab_equation,
      hjust = 0.5, vjust = 0.5,
      size = 4
    ) +
    labs(
      x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      y = "PTZ Width (km)"
    ) +
    theme_bw(base_size = 14) +
    theme_2() +
    theme(
      plot.tag.position = "topright",
      plot.tag = element_text(margin = margin(5, 10, 0, 0))
    )

  p4 <-
    df_plume |>
    ggplot(aes(x = max_reaction_rate, y = displacement)) +
    geom_line(
      data = log_plume_predictions2, aes(x = max_reaction_rate, y = displacement),
      color = "black", linewidth = 0.8, inherit.aes = FALSE
    ) +
    geom_point(color = "black") +
    annotate("text",
      x = ((max(df_plume$max_reaction_rate) - min(df_plume$max_reaction_rate)) / 2) + 2,
      y = ((max(df_plume$displacement) - min(df_plume$displacement)) / 2) - 13,
      label = log_plume_equation2,
      hjust = 0.5, vjust = 0.5,
      size = 4
    ) +
    labs(
      x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      y = "PTZ Displacement (km)"
    ) +
    theme_bw(base_size = 14) +
    theme_2() +
    theme(
      plot.tag.position = "topright",
      plot.tag = element_text(margin = margin(5, 10, 0, 0))
    )

  p5 <-
    filter(df_slab, displacement > -85) |>
    ggplot(aes(x = max_reaction_rate, y = displacement)) +
    geom_line(
      data = log_slab_predictions, aes(x = max_reaction_rate, y = displacement),
      color = "black", linetype = 2, linewidth = 0.5, inherit.aes = FALSE
    ) +
    geom_line(
      data = log_slab_predictions_no_anomalies, aes(x = max_reaction_rate, y = displacement),
      color = "black", linewidth = 0.8, inherit.aes = FALSE
    ) +
    geom_point(data = filter(df_slab, displacement <= -85), color = "red") +
    geom_point(color = "black") +
    annotate("text",
      x = 4.4, y = -100,
      label = log_slab_equation,
      hjust = 1, vjust = 0,
      size = 4
    ) +
    scale_y_continuous(limits = c(NA, 47)) +
    labs(
      x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      y = "PTZ Displacement (km)",
    ) +
    theme_bw(base_size = 14) +
    theme_2() +
    theme(
      plot.tag.position = "topright",
      plot.tag = element_text(margin = margin(5, 10, 0, 0))
    )

  p <- p0 + p2 + p4 + plot_annotation(tag_levels = "a")

  ggsave(
    out_path1,
    plot = p,
    width = 8.6,
    height = 3.0,
    dpi = 300,
    bg = "white"
  )

  p <- p1 + p3 + p5 + plot_annotation(tag_levels = "a")

  ggsave(
    out_path2,
    plot = p,
    width = 8.6,
    height = 3.0,
    dpi = 300,
    bg = "white"
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript main.R [in_dir] [out_dir]\n")
    cat(" !! Example: Rscript main.R /path/to/sim_out_dir /path/to/fig_dir\n")
    return(invisible(NULL))
  }

  in_dir <- args[1]
  out_dir <- args[2]

  in_displacement <- file.path(in_dir, "centerline-profile-data.csv")
  out_plumes <- file.path(out_dir, "ptz-plumes.png")
  out_slabs <- file.path(out_dir, "ptz-slabs.png")
  visualize_displacement(in_displacement, out_plumes, out_slabs)
}

if (
  !interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))
) {
  main()
}
