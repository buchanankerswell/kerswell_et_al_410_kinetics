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
## .1. Visualize Statistics                      !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_statistics <- function(statistics_path, residuals_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  stats <- read_aspect_stats(statistics_path) |> filter(`time_(years)` > 5e6)
  res <- read_stokes_residuals(residuals_path) |>
    filter(time > 5e6) |>
    mutate(
      running_mean = rollmean(mean_residual, 20, fill = NA, align = "center")
    )

  p0 <-
    ggplot(mapping = aes(x = time / 1e6)) +
    geom_path(data = res, aes(y = mean_residual), alpha = 0.1) +
    geom_path(data = res, aes(y = running_mean), linewidth = 1.0) +
    labs(
      x = "Time (Ma)",
      y = NULL,
      title = "Avgerage Stokes residuals (log10)"
    ) +
    theme_bw(base_size = 14) +
    theme_1() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  p1 <-
    ggplot(mapping = aes(
      x = `time_(years)` / 1e6,
      y = `rms_velocity_(m/year)` * 1e2
    )) +
    geom_path(data = stats, linewidth = 1.0) +
    labs(
      x = "Time (Ma)",
      y = NULL,
      title = expression(paste("RMS ", italic(u), " (cm/yr)"))
    ) +
    theme_bw(base_size = 14) +
    theme_1() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  p2 <-
    ggplot(mapping = aes(
      x = `time_(years)` / 1e6,
      y = `max._velocity_(m/year)` * 1e2
    )) +
    geom_path(data = stats, linewidth = 1.0) +
    labs(
      x = "Time (Ma)",
      y = NULL,
      title = expression(paste("Max ", italic(u), " (cm/yr)"))
    ) +
    theme_bw(base_size = 14) +
    theme_1() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  p3 <-
    ggplot(mapping = aes(
      x = `time_(years)` / 1e6,
      y = `average_temperature_(k)`
    )) +
    geom_path(data = stats, linewidth = 1.0) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    labs(x = "Time (Ma)", y = NULL, title = "Average T (K)") +
    theme_bw(base_size = 14) +
    theme_1()

  p <- p0 / p1 / p2 / p3

  ggsave(
    out_path,
    plot = p,
    width = 6.5,
    height = 8.5,
    dpi = 300,
    bg = "white"
  )
}
