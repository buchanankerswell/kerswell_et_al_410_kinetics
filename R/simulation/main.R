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
visualize_displacement <- function(in_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  if (!file.exists(in_path)) {
    cat(" !! Missing data: ", basename(in_path), "\n", sep = "")
    return(invisible())
  }

  df <- read_displacements(in_path) |>
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

  width_range <- range(df_slab$width, df_plume$width)
  displacement_range <- range(df_slab$displacement, df_plume$displacement)

  p0 <- df_plume |>
    ggplot(aes(x = displacement, y = width, fill = max_reaction_rate)) +
    geom_point(size = 3.0, shape = 21, color = "black") +
    scale_y_continuous(limits = width_range, expand = expansion(mult = c(0.1, 0.1))) +
    scale_fill_viridis_c(
      name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      trans = "log10", option = "plasma",
      breaks = c(1e-1, 1e1, 1e3, 1e5),
      labels = label_log()
    ) +
    labs(x = "Displacement (km)", y = "Width (km)", title = "410 Plumes") +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(plot.tag.position = "topright", plot.tag = element_text(margin = margin(5, 10, 0, 0)))

  p1 <- df_slab |>
    ggplot(aes(x = displacement, y = width, fill = max_reaction_rate)) +
    geom_point(size = 3.0, shape = 21, color = "black") +
    scale_y_continuous(limits = width_range, expand = expansion(mult = c(0.1, 0.1))) +
    scale_fill_viridis_c(
      name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      trans = "log10", option = "plasma",
      breaks = c(1e-1, 1e1, 1e3),
      labels = label_log()
    ) +
    labs(x = "Displacement (km)", y = NULL, title = "410 Slabs") +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(plot.tag.position = "topright", plot.tag = element_text(margin = margin(5, 10, 0, 0)))

  p2 <- df_plume |>
    ggplot(aes(x = max_reaction_rate, y = width, fill = max_reaction_rate)) +
    geom_point(size = 3.0, shape = 21, color = "black", show.legend = FALSE) +
    scale_x_continuous(trans = "log10", labels = label_log()) +
    scale_y_continuous(limits = width_range, expand = expansion(mult = c(0.1, 0.1))) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      trans = "log10", option = "plasma",
      breaks = c(1e-1, 1e1, 1e3, 1e5),
      labels = label_log()
    ) +
    labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = "Width (km)") +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(plot.tag.position = "topright", plot.tag = element_text(margin = margin(5, 10, 0, 0)))

  p3 <- df_slab |>
    ggplot(aes(x = max_reaction_rate, y = width, fill = max_reaction_rate)) +
    geom_point(size = 3.0, shape = 21, color = "black", show.legend = FALSE) +
    scale_x_continuous(trans = "log10", labels = label_log()) +
    scale_y_continuous(limits = width_range, expand = expansion(mult = c(0.1, 0.1))) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      trans = "log10", option = "plasma",
      breaks = c(1e-1, 1e1, 1e3),
      labels = label_log()
    ) +
    labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = NULL) +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(plot.tag.position = "topright", plot.tag = element_text(margin = margin(5, 10, 0, 0)))

  p4 <- df_plume |>
    ggplot(aes(x = max_reaction_rate, y = displacement, fill = max_reaction_rate)) +
    geom_point(size = 3.0, shape = 21, color = "black", show.legend = FALSE) +
    scale_x_continuous(trans = "log10", labels = label_log()) +
    scale_y_continuous(limits = displacement_range, expand = expansion(mult = c(0.1, 0.1))) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      trans = "log10", option = "plasma",
      breaks = c(1e-1, 1e1, 1e3, 1e5),
      labels = label_log()
    ) +
    labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = "Displacement (km)") +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(plot.tag.position = "topright", plot.tag = element_text(margin = margin(5, 10, 0, 0)))

  p5 <- df_slab |>
    ggplot(aes(x = max_reaction_rate, y = displacement, fill = max_reaction_rate)) +
    geom_point(size = 3.0, shape = 21, color = "black", show.legend = FALSE) +
    scale_x_continuous(trans = "log10", labels = label_log()) +
    scale_y_continuous(limits = displacement_range, expand = expansion(mult = c(0.1, 0.1))) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
      trans = "log10", option = "plasma",
      breaks = c(1e-1, 1e1, 1e3),
      labels = label_log()
    ) +
    labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = NULL, ) +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(plot.tag.position = "topleft", plot.tag = element_text(margin = margin(5, 0, 0, 0)))

  pp1 <- p0 / p2 / p4 + plot_annotation(tag_levels = "a") & theme(plot.title = element_text(size = 16, hjust = 0.5))
  pp2 <- p1 / p3 / p5 + plot_annotation(tag_levels = list(c("d", "e", "f"))) & theme(plot.title = element_text(size = 16, hjust = 0.5))
  p <- pp1 | pp2

  ggsave(out_path, plot = p, width = 6.0, height = 8.6, dpi = 300, bg = "white")
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

  in_displacement <- file.path(in_dir, "depth-profile-data.csv")
  out_comp <- file.path(out_dir, "410-structure.png")

  cat("    --------------------------------------------------\n")
  cat("    Drawing depth profile summary\n")
  cat("    --------------------------------------------------\n")

  visualize_displacement(in_displacement, out_comp)
}

if (
  !interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))
) {
  main()
}
