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
## .1. Visualize 410 structure                   !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_displacement <- function(in_path, out_path) {
  if (!file.exists(in_path)) {
    cat(" !! Missing data: ", basename(in_path), "\n", sep = "")
    return(invisible())
  }

  df <- read_displacements(in_path) |>
    filter(timestep == 10) |>
    mutate(model_id = str_replace_all(model_id, "_", "-"))
  df_slab <- df |>
    filter(str_detect(model_id, "slab")) |>
    mutate(model_id = str_replace_all(model_id, "slab-", "")) |>
    mutate(model_id = str_replace_all(model_id, "-", " "))
  df_plume <- df |>
    filter(str_detect(model_id, "plume")) |>
    mutate(model_id = str_replace_all(model_id, "plume-", "")) |>
    mutate(model_id = str_replace_all(model_id, "-", " ")) |>
    mutate(width = abs(width))

  calculate_tile_boundaries <- function(data_frame) {
    z_unique <- sort(unique(data_frame$Z_factor))

    z_boundaries <- c(
      z_unique[1] / sqrt(z_unique[2] / z_unique[1]),
      sqrt(z_unique[-length(z_unique)] * z_unique[-1]),
      z_unique[length(z_unique)] * sqrt(z_unique[length(z_unique)] / z_unique[length(z_unique) - 1])
    )

    z_factor_map <- data.frame(
      Z_factor = z_unique,
      xmin = z_boundaries[-length(z_boundaries)],
      xmax = z_boundaries[-1]
    )

    data_frame |> left_join(z_factor_map, by = "Z_factor")
  }

  df_slab_tiled <- calculate_tile_boundaries(df_slab)
  df_plume_tiled <- calculate_tile_boundaries(df_plume)

  width_range <- range(df_slab$width, df_plume$width)
  displacement_range <- range(df_slab$displacement, df_plume$displacement)
  max_velocity_range <- range(df_slab$max_velocity, df_plume$max_velocity)
  rate_range <- range(df_slab$max_reaction_rate, df_plume$max_reaction_rate)
  z_range <- range(df_slab_tiled$xmin, df_slab_tiled$xmax, df_plume_tiled$xmin, df_plume_tiled$xmax)

  plume_breaks_rxn <- c(1e-1, 1e1, 1e3, 1e5)
  slab_breaks_rxn <- c(1e-1, 1e1, 1e3, 1e5)
  breaks_z <- c(1e1, 1e3, 1e5, 1e7)

  draw_composition <- function(b) {
    out_pth <- str_replace(out_path, ".png", paste0("-B", b, ".png"))

    if (plot_exists(out_pth)) {
      return(invisible())
    }

    d_slab <- df_slab |> filter(B_factor == b)
    d_plume <- df_plume |> filter(B_factor == b)

    p0 <- d_plume |>
      ggplot(aes(x = max_reaction_rate, y = width, fill = max_reaction_rate)) +
      geom_point(size = 4, color = "black", shape = 21) +
      scale_x_continuous(trans = "log10", labels = label_log(), breaks = plume_breaks_rxn, expand = expansion(mult = c(0.1, 0.1))) +
      scale_y_continuous(limits = width_range, expand = expansion(mult = c(0.1, 0.1))) +
      annotation_logticks(sides = "b", linewidth = 0.2) +
      scale_fill_viridis_c(
        name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
        trans = "log10", option = "plasma",
        breaks = plume_breaks_rxn,
        labels = label_log(),
        limits = rate_range
      ) +
      labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = "Width (km)", title = paste0("Plumes: B = ", b)) +
      theme_bw(base_size = 16) +
      theme_2() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())

    p1 <- d_slab |>
      ggplot(aes(x = max_reaction_rate, y = width, fill = max_reaction_rate)) +
      geom_point(size = 4, color = "black", shape = 21, show.legend = FALSE) +
      scale_x_continuous(trans = "log10", labels = label_log(), breaks = slab_breaks_rxn, expand = expansion(mult = c(0.1, 0.1))) +
      scale_y_continuous(limits = width_range, expand = expansion(mult = c(0.1, 0.1))) +
      annotation_logticks(sides = "b", linewidth = 0.2) +
      scale_fill_viridis_c(
        name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
        trans = "log10", option = "plasma",
        breaks = slab_breaks_rxn,
        labels = label_log(),
        limits = rate_range
      ) +
      labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = NULL, title = paste0("Slabs: B = ", b)) +
      theme_bw(base_size = 16) +
      theme_2() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

    p2 <- d_plume |>
      ggplot(aes(x = max_reaction_rate, y = displacement, fill = max_reaction_rate)) +
      geom_point(size = 4, color = "black", shape = 21, show.legend = FALSE) +
      scale_x_continuous(trans = "log10", labels = label_log(), breaks = plume_breaks_rxn, expand = expansion(mult = c(0.1, 0.1))) +
      scale_y_reverse(limits = displacement_range, expand = expansion(mult = c(0.1, 0.1))) +
      annotation_logticks(sides = "b", linewidth = 0.2) +
      scale_fill_viridis_c(
        name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
        trans = "log10", option = "plasma",
        breaks = plume_breaks_rxn,
        labels = label_log(),
        limits = rate_range
      ) +
      labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = "Displacement (km)") +
      theme_bw(base_size = 16) +
      theme_2()

    p3 <- d_slab |>
      ggplot(aes(x = max_reaction_rate, y = displacement, fill = max_reaction_rate)) +
      geom_point(size = 4, color = "black", shape = 21, show.legend = FALSE) +
      scale_x_continuous(trans = "log10", labels = label_log(), breaks = slab_breaks_rxn, expand = expansion(mult = c(0.1, 0.1))) +
      scale_y_reverse(limits = displacement_range, expand = expansion(mult = c(0.1, 0.1))) +
      annotation_logticks(sides = "b", linewidth = 0.2) +
      scale_fill_viridis_c(
        name = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"),
        trans = "log10", option = "plasma",
        breaks = slab_breaks_rxn,
        labels = label_log(),
        limits = rate_range
      ) +
      labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = NULL) +
      theme_bw(base_size = 16) +
      theme_2() +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank())

    p <- ((p0 / p2) | (p1 / p3)) & theme(plot.title = element_text(size = 16))
    ggsave(out_pth, plot = p, width = 5.0, height = 5.0, dpi = 300, bg = "white")
  }

  walk(unique(c(df_slab$B_factor, df_plume$B_factor)), draw_composition)

  draw_rect <- function(b_row, zmin, zmax, color = "black", size = 0.3, linetype = 1, alpha = 1.0) {
    z_unique <- sort(unique(df_slab$Z_factor))
    z_step_factor <- sqrt(z_unique[2] / z_unique[1])
    geom_rect(
      aes(xmin = zmin / z_step_factor, xmax = zmax * z_step_factor, ymin = b_row - 1.0, ymax = b_row + 1.0),
      fill = NA,
      color = color,
      alpha = alpha,
      linewidth = size,
      linetype = linetype
    )
  }

  rects_slab_b <- list(
    draw_rect(2, 1.8e5, 7.0e7),
    draw_rect(2, 2.0e2, 7.8e4),
    draw_rect(2, 3.0e0, 8.7e1),
    draw_rect(4, 7.8e4, 7.0e7),
    draw_rect(4, 2.0e2, 3.3e4),
    draw_rect(4, 3.0e0, 8.7e1),
    draw_rect(6, 3.3e4, 7.0e7),
    draw_rect(6, 3.7e1, 1.4e4),
    draw_rect(6, 3.0e0, 1.6e1),
    draw_rect(8, 6.0e3, 7.0e7),
    draw_rect(8, 1.6e1, 2.6e3),
    draw_rect(8, 3.0e0, 7.0e0),
    draw_rect(10, 6.0e3, 7.0e7),
    draw_rect(10, 1.6e1, 2.6e3),
    draw_rect(10, 3.0e0, 7.0e0)
  )
  rects_slab_w <- list(
    draw_rect(2, 1.8e5, 7.0e7, color = "white"),
    draw_rect(2, 2.0e2, 7.8e4, color = "white"),
    draw_rect(2, 3.0e0, 8.7e1, color = "white"),
    draw_rect(4, 7.8e4, 7.0e7, color = "white"),
    draw_rect(4, 2.0e2, 3.3e4, color = "white"),
    draw_rect(4, 3.0e0, 8.7e1, color = "white"),
    draw_rect(6, 3.3e4, 7.0e7, color = "white"),
    draw_rect(6, 3.7e1, 1.4e4, color = "white"),
    draw_rect(6, 3.0e0, 1.6e1, color = "white"),
    draw_rect(8, 6.0e3, 7.0e7, color = "white"),
    draw_rect(8, 1.6e1, 2.6e3, color = "white"),
    draw_rect(8, 3.0e0, 7.0e7, color = "white"),
    draw_rect(10, 6.0e3, 7.0e7, color = "white"),
    draw_rect(10, 1.6e1, 2.6e3, color = "white"),
    draw_rect(10, 3.0e0, 7.0e0, color = "white")
  )
  rects_plume_b <- list(
    draw_rect(2, 3.0e0, 7.0e7),
    draw_rect(4, 3.0e0, 7.0e7),
    draw_rect(6, 3.0e0, 7.0e7),
    draw_rect(8, 3.0e0, 7.0e7),
    draw_rect(10, 3.0e0, 7.0e7)
  )
  rects_plume_w <- list(
    draw_rect(2, 3.0e0, 7.0e7, color = "white"),
    draw_rect(4, 3.0e0, 7.0e7, color = "white"),
    draw_rect(6, 3.0e0, 7.0e7, color = "white"),
    draw_rect(8, 3.0e0, 7.0e7, color = "white"),
    draw_rect(10, 3.0e0, 7.0e7, color = "white")
  )

  p0 <- df_slab_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = displacement)) +
    geom_rect(show.legend = FALSE) +
    rects_slab_b +
    annotate("text", x = 2.4e6, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "equilibrium", fontface = "bold") +
    annotate("text", x = 2.6e3, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "int.", fontface = "bold") +
    annotate("text", x = 1.6e1, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "slug.", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_gradient2(low = "#A50026", mid = "white", high = "#313695", limits = displacement_range) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = "Displacement (km)") +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())

  p1 <- df_slab_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = width)) +
    geom_rect(show.legend = FALSE) +
    rects_slab_b +
    annotate("text", x = 2.4e6, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "sharp", fontface = "bold") +
    annotate("text", x = 2.6e3, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "fading", fontface = "bold") +
    annotate("text", x = 1.6e1, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "sharp", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(option = "mako", limits = width_range, direction = -1) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = "Width (km)", title = "Slab Regimes") +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())

  p2 <- df_slab_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = max_velocity)) +
    geom_rect(show.legend = FALSE) +
    rects_slab_w +
    annotate("text", x = 2.4e6, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "penetrate", color = "white", fontface = "bold") +
    annotate("text", x = 2.6e3, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "int.", color = "white", fontface = "bold") +
    annotate("text", x = 1.6e1, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "pond", color = "white", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(option = "rocket", limits = max_velocity_range) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = bquote("Max " * italic(u)[y] * " (cm/yr)")) +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())

  p3 <- df_plume_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = displacement)) +
    geom_rect() +
    rects_plume_b +
    annotate("text", x = 1.4e4, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "equilibrium", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_gradient2(low = "#A50026", mid = "white", high = "#313695", limits = displacement_range) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = "Displacement (km)") +
    theme_bw(base_size = 16) +
    theme_2()

  p4 <- df_plume_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = width)) +
    geom_rect() +
    rects_plume_b +
    annotate("text", x = 1.4e4, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "sharp", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(option = "mako", limits = width_range, direction = -1) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = "Width (km)", title = "Plume Regimes") +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank())

  p5 <- df_plume_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = max_velocity)) +
    geom_rect() +
    rects_plume_b +
    annotate("text", x = 1.4e4, y = 4, size = 3, hjust = 0.5, vjust = 0.5, label = "penetrate", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(option = "rocket", limits = max_velocity_range) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = bquote("Max " * italic(u)[y] * " (cm/yr)")) +
    theme_bw(base_size = 16) +
    theme_2() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.title = element_text(margin = margin(0, 0, 2, 0)))

  p <- (p0 | p1 | p2) / (p3 | p4 | p5) & theme(plot.title = element_text(size = 16))

  out_path_comp <- str_replace(out_path, ".png", "-comp.png")

  if (plot_exists(out_path_comp)) {
    return(invisible())
  }

  suppressWarnings(ggsave(out_path_comp, plot = p, width = 6.5, height = 5.0, dpi = 300, bg = "white"))

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

  in_data <- file.path(in_dir, "depth-profile-data.csv")
  out_path <- file.path(out_dir, "410-structure.png")

  cat("    --------------------------------------------------\n")
  cat("    Drawing 410 structure summary\n")
  cat("    --------------------------------------------------\n")

  visualize_displacement(in_data, out_path)
}

if (
  !interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))
) {
  main()
}
