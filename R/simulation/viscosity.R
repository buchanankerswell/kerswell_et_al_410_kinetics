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
## .1. Visualize Viscosity Profiles              !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_viscosity_profile <- function(profile_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  profile <- read_depth_average(profile_path)

  middle_time <-
    profile$time[which.min(abs(profile$time - max(profile$time) / 2))]

  p0 <-
    profile |>
    filter(time == 0) |>
    ggplot(mapping = aes(log_viscosity, depth / 1e3)) +
    geom_path(linewidth = 1.0) +
    scale_x_continuous(limits = c(19, 24)) +
    scale_y_reverse() +
    labs(x = bquote("Log" ~ eta ~ (Pa ~ s)), y = "Depth (km)") +
    annotate("text",
      x = 19, y = Inf, label = bquote(eta[bm](t[0])),
      vjust = -1.5, hjust = 0, size = 7
    ) +
    theme_bw(base_size = 22) +
    theme_1()

  p1 <-
    profile |>
    filter(time == middle_time) |>
    ggplot(mapping = aes(log_viscosity, depth / 1e3)) +
    geom_path(linewidth = 1.0) +
    scale_x_continuous(limits = c(19, 24)) +
    scale_y_reverse() +
    labs(x = bquote("Log" ~ eta ~ (Pa ~ s)), y = "Depth (km)") +
    annotate("text",
      x = 19, y = Inf, label = bquote(eta[bm](t[int])),
      vjust = -1.5, hjust = 0, size = 7
    ) +
    theme_bw(base_size = 22) +
    theme_1() +
    theme(axis.title.y = element_blank())

  p2 <-
    profile |>
    filter(time == max(time)) |>
    ggplot(mapping = aes(log_viscosity, depth / 1e3)) +
    geom_path(linewidth = 1.0) +
    scale_x_continuous(limits = c(19, 24)) +
    scale_y_reverse() +
    labs(x = bquote("Log" ~ eta ~ (Pa ~ s)), y = "Depth (km)") +
    annotate("text",
      x = 19, y = Inf, label = bquote(eta[bm](t[end])),
      vjust = -1.5, hjust = 0, size = 7
    ) +
    theme_bw(base_size = 22) +
    theme_1() +
    theme(axis.title.y = element_blank())

  p <- (p0 | p1 | p2) + plot_annotation(tag_levels = "a")

  ggsave(
    out_path,
    plot = p,
    width = 10.5,
    height = 5.5,
    dpi = 300,
    bg = "white"
  )
}
