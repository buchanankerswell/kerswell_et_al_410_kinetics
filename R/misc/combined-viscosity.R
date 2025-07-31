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
visualize_viscosity_profiles <- function(profile_paths, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  parent_dirs <- map_chr(profile_paths, ~ basename(dirname(.x)))
  profile_names <- map_chr(parent_dirs, ~ {
    str_replace_all(.x, "_", " ")
  })

  all_times <-
    map(profile_paths, ~ {
      read_depth_average(.x)$time
    }) |>
    unlist()

  min_time <- min(all_times, na.rm = TRUE)
  max_time <- max(all_times, na.rm = TRUE)

  viscosity_init <-
    map2(profile_paths, profile_names, ~ {
      read_depth_average(.x) |>
        select(c(time, depth, log_viscosity)) |>
        filter(time == 0) |>
        select(-time) |>
        rename(!!paste0("log_viscosity_", .y) := log_viscosity)
    }) |>
    reduce(full_join, by = "depth") |>
    pivot_longer(
      cols = starts_with("log_viscosity_"),
      names_to = "name",
      values_to = "log_viscosity"
    ) |>
    mutate(name = gsub("log_viscosity_", "", name)) |>
    mutate(name = factor(name, levels = profile_names))

  viscosity_end <-
    map2(profile_paths, profile_names, ~ {
      read_depth_average(.x) |>
        select(c(time, depth, log_viscosity)) |>
        filter(time == max(time)) |>
        select(-time) |>
        rename(!!paste0("log_viscosity_", .y) := log_viscosity)
    }) |>
    reduce(full_join, by = "depth") |>
    pivot_longer(
      cols = starts_with("log_viscosity_"),
      names_to = "name",
      values_to = "log_viscosity"
    ) |>
    mutate(name = gsub("log_viscosity_", "", name)) |>
    mutate(name = factor(name, levels = profile_names))

  p0 <-
    ggplot(
      viscosity_init,
      aes(x = log_viscosity, y = depth / 1e6, color = name)
    ) +
    geom_path(linewidth = 0.7) +
    annotate(
      "text",
      x = -Inf,
      y = 2,
      label = paste0(round(min_time, 1), " Ma"),
      hjust = -0.3,
      vjust = 0,
      size = 5
    ) +
    scale_y_reverse(expand = expansion(mult = 0.15)) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = bquote("Log" ~ eta ~ "[Pa s]"),
      y = bquote("Depth [" * km %*% 10^3 * "]"),
      color = NULL
    ) +
    theme_bw(base_size = 18) +
    profile_theme()

  p1 <-
    ggplot(
      viscosity_end,
      aes(x = log_viscosity, y = depth / 1e6, color = name)
    ) +
    geom_path(linewidth = 0.7, show.legend = FALSE) +
    annotate(
      "text",
      x = -Inf,
      y = 2,
      label = paste0(round(max_time, 1), " Ma"),
      hjust = -0.3,
      vjust = 0,
      size = 5
    ) +
    scale_y_reverse(expand = expansion(mult = 0.15)) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = bquote("Log" ~ eta ~ "[Pa s]"),
      y = bquote("Depth [" * km %*% 10^3 * "]"),
      color = NULL
    ) +
    theme_bw(base_size = 18) +
    profile_theme() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    )

  p <- (p0 + p1) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom")

  ggsave(
    out_path,
    plot = p,
    width = 6.5,
    height = 4.5,
    dpi = 300,
    bg = "white"
  )
}
