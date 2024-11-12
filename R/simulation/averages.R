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
## .1. Visualize Depth Averages                  !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_depth_averages <- function(profile_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  profile <-
    read_depth_average(profile_path) |>
    mutate(sinking_velocity = sinking_velocity * 10)

  variable_order <- c(
    "temperature",
    "log_viscosity",
    "sinking_velocity",
    "vertical_heat_flux"
  )

  df_long <-
    profile |>
    pivot_longer(
      cols = all_of(variable_order),
      names_to = "variable",
      values_to = "value"
    )

  df_long$variable <- factor(df_long$variable, levels = variable_order)

  units <- c(
    "temperature" = "K",
    "log_viscosity" = "log~Pa~s",
    "sinking_velocity" = "mm/yr",
    "vertical_heat_flux" = "mW/m^2"
  )

  labs <- c(
    "temperature" = "T",
    "log_viscosity" = "eta",
    "sinking_velocity" = "u[sink]",
    "vertical_heat_flux" = "Q"
  )

  custom_labeller <- function(variable) {
    paste0(labs[variable], " * ' [' * ", units[variable], " * ']'")
  }

  p <-
    df_long |>
    ggplot(aes(x = value, y = depth / 1e6)) +
    geom_path(
      data = filter(df_long, time == max(time)),
      aes(linetype = paste0(round(max(time), 1), " Ma")),
      linewidth = 1.5
    ) +
    geom_path(
      data = filter(df_long, time == min(time)),
      aes(linetype = "0 Ma"),
      linewidth = 1.5
    ) +
    scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.15)) +
    scale_y_reverse(expand = expansion(mult = 0.15)) +
    labs(x = NULL, y = bquote("Depth [" * km %*% 10^3 * "]"), linetype = NULL) +
    facet_wrap(~variable,
      scales = "free_x", nrow = 1,
      labeller = labeller(
        variable = custom_labeller,
        .default = label_parsed
      )
    ) +
    theme_bw(base_size = 36) +
    profile_theme() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold"),
      legend.direction = "vertical",
      legend.position.inside = c(0.01, 0.15)
    )

  ggsave(
    out_path,
    plot = p,
    width = 20,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}
