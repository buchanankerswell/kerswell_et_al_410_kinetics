#######################################################
## Visualize 410 structure                       !!! ##
#######################################################
# ASPECT simulations ran with "Use years instead of seconds = true"
# Stored (per-year) Z values are converted to per-second by dividing by this factor
SECONDS_PER_YEAR <- 3.15e7

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_410_data <- function(filepath) {
  suppressWarnings(df <- read_csv(filepath, show_col_types = FALSE))
  df |>
    mutate(Z_factor = Z_factor / SECONDS_PER_YEAR) |>
    mutate(displacement = displacement / 1e3, width = width / 1e3)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_410_structure <- function(in_path, out_path) {
  if (!file.exists(in_path)) {
    cat(" !! Missing data: ", basename(in_path), "\n", sep = "")
    return(invisible())
  }

  df <- read_410_data(in_path) |>
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

  width_range_slab <- range(df_slab$width)
  displacement_range_slab <- range(df_slab$displacement)
  max_velocity_range_slab <- range(df_slab$max_velocity)
  rate_range_slab <- range(df_slab$max_reaction_rate)
  z_range_slab <- range(df_slab_tiled$xmin, df_slab_tiled$xmax)

  plume_breaks_rxn <- c(1e-1, 1e1, 1e3, 1e5)
  slab_breaks_rxn <- c(1e-1, 1e1, 1e3, 1e5)
  breaks_z <- c(1e-7, 1e-5, 1e-3, 1e-1)

  draw_composition <- function(b) {
    out_pth <- str_replace(out_path, ".png", paste0("-B", b, ".png"))

    if (plot_exists(out_pth)) {
      return(invisible())
    }

    d_slab <- df_slab |> filter(B_factor == b)
    d_plume <- df_plume |> filter(B_factor == b)

    p0 <- d_plume |>
      ggplot(aes(x = max_reaction_rate, y = width, fill = max_reaction_rate)) +
      geom_point(size = 3.5, color = "black", shape = 21) +
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
      labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = "Width (km)", title = bquote("Plumes: " * italic(B) * " = " * .(b))) +
      theme_bw(base_size = 14) +
      theme_1() +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(margin = margin(0, 0, 2, 0)),
        legend.margin = margin(-6, 0, 0, 0)
      )

    p1 <- d_slab |>
      ggplot(aes(x = max_reaction_rate, y = width, fill = max_reaction_rate)) +
      geom_hline(yintercept = 173, linewidth = 0.3, color = "grey50") +
      geom_vline(xintercept = 2.14e0, linewidth = 0.3, color = "grey50") +
      geom_vline(xintercept = 6.76e-2, linewidth = 0.3, color = "grey50") +
      geom_point(size = 3.5, color = "black", shape = 21, show.legend = FALSE) +
      annotate("text", x = 2.3e1, y = 180, label = "(1)", size = 4, vjust = 0, fontface = "bold") +
      annotate("text", x = 3.8e-1, y = 180, label = "(2)", size = 4, vjust = 0, fontface = "bold") +
      annotate("text", x = 1.8e-2, y = 180, label = "(3)", size = 4, vjust = 0, fontface = "bold") +
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
      labs(x = bquote(dot(italic(X))["max"] * " (" * Ma^-1 * ")"), y = NULL, title = bquote("Slabs: " * italic(B) * " = " * .(b))) +
      theme_bw(base_size = 14) +
      theme_1() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

    p2 <- d_plume |>
      ggplot(aes(x = max_reaction_rate, y = displacement, fill = max_reaction_rate)) +
      geom_point(size = 3.5, color = "black", shape = 21, show.legend = FALSE) +
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
      theme_bw(base_size = 14) +
      theme_1()

    p3 <- d_slab |>
      ggplot(aes(x = max_reaction_rate, y = displacement, fill = max_reaction_rate)) +
      geom_vline(xintercept = 2.14e0, linewidth = 0.3, color = "grey50") +
      geom_vline(xintercept = 6.76e-2, linewidth = 0.3, color = "grey50") +
      geom_point(size = 3.5, color = "black", shape = 21, show.legend = FALSE) +
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
      theme_bw(base_size = 14) +
      theme_1() +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank())

    p <- (p0 / p2) | (p1 / p3)
    ggsave(out_pth, plot = p, width = 4.5, height = 4.2, dpi = 300, bg = "white")
  }

  draw_composition(b = 4)

  draw_rect <- function(b_row, zmin, zmax, color = "black", size = 0.3, linetype = 1, alpha = 1.0) {
    z_unique <- sort(unique(df_slab$Z_factor))
    z_step_factor <- sqrt(z_unique[2] / z_unique[1])

    # Clamp values to dataset boundaries
    xmin_val <- max(zmin / z_step_factor, min(df_slab_tiled$xmin))
    xmax_val <- min(zmax * z_step_factor, max(df_slab_tiled$xmax))

    geom_rect(
      aes(xmin = xmin_val, xmax = xmax_val, ymin = b_row - 1.0, ymax = b_row + 1.0),
      fill = NA,
      color = color,
      alpha = alpha,
      linewidth = size,
      linetype = linetype
    )
  }

  rects_slab_b <- list(
    draw_rect(2, 1.4e-2, 2.2e0),
    draw_rect(2, 6.3e-6, 5.7e-3),
    draw_rect(2, 9.5e-8, 2.8e-6),
    draw_rect(4, 5.7e-3, 2.2e0),
    draw_rect(4, 6.3e-6, 2.5e-3),
    draw_rect(4, 9.5e-8, 2.8e-6),
    draw_rect(6, 2.5e-3, 2.2e0),
    draw_rect(6, 1.2e-6, 1.0e-3),
    draw_rect(6, 9.5e-8, 5.1e-7),
    draw_rect(8, 4.4e-4, 2.2e0),
    draw_rect(8, 5.1e-7, 1.9e-4),
    draw_rect(8, 9.5e-8, 2.2e-7),
    draw_rect(10, 4.4e-4, 2.2e0),
    draw_rect(10, 5.1e-7, 1.9e-4),
    draw_rect(10, 9.5e-8, 2.2e-7)
  )
  rects_slab_w <- list(
    draw_rect(2, 1.4e-2, 2.2e0, color = "white"),
    draw_rect(2, 6.3e-6, 5.7e-3, color = "white"),
    draw_rect(2, 9.5e-8, 2.8e-6, color = "white"),
    draw_rect(4, 5.7e-3, 2.2e0, color = "white"),
    draw_rect(4, 6.3e-6, 2.5e-3, color = "white"),
    draw_rect(4, 9.5e-8, 2.8e-6, color = "white"),
    draw_rect(6, 2.5e-3, 2.2e0, color = "white"),
    draw_rect(6, 1.2e-6, 1.0e-3, color = "white"),
    draw_rect(6, 9.5e-8, 5.1e-7, color = "white"),
    draw_rect(8, 4.4e-4, 2.2e0, color = "white"),
    draw_rect(8, 5.1e-7, 1.9e-4, color = "white"),
    draw_rect(8, 9.5e-8, 2.2e-7, color = "white"),
    draw_rect(10, 4.4e-4, 2.2e0, color = "white"),
    draw_rect(10, 5.1e-7, 1.9e-4, color = "white"),
    draw_rect(10, 9.5e-8, 2.2e-7, color = "white")
  )
  rects_plume_b <- list(
    draw_rect(2, 9.5e-8, 2.2e0),
    draw_rect(4, 9.5e-8, 2.2e0),
    draw_rect(6, 9.5e-8, 2.2e0),
    draw_rect(8, 9.5e-8, 2.2e0),
    draw_rect(10, 9.5e-8, 2.2e0)
  )
  rects_plume_w <- list(
    draw_rect(2, 9.5e-8, 2.2e0, color = "white"),
    draw_rect(4, 9.5e-8, 2.2e0, color = "white"),
    draw_rect(6, 9.5e-8, 2.2e0, color = "white"),
    draw_rect(8, 9.5e-8, 2.2e0, color = "white"),
    draw_rect(10, 9.5e-8, 2.2e0, color = "white")
  )

  p0 <- df_slab_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = max_reaction_rate)) +
    geom_rect() +
    rects_slab_w +
    annotate("text", x = 1.1e-1, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(1)", fontface = "bold", color = "white") +
    annotate("text", x = 1.3e-4, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(2)", fontface = "bold", color = "white") +
    annotate("text", x = 5.1e-7, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(3)", fontface = "bold", color = "white") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range_slab, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(trans = "log10", option = "viridis", labels = label_log(), limits = rate_range_slab) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = bquote("Max " * italic(dot(X)) * " (Ma"^-1 * ")")) +
    theme_bw(base_size = 14) +
    theme_2() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.box.margin = margin(0, 0, 2, 0))

  p1 <- df_slab_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = displacement)) +
    geom_rect() +
    rects_slab_b +
    annotate("text", x = 1.1e-1, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(1)", fontface = "bold") +
    annotate("text", x = 1.3e-4, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(2)", fontface = "bold") +
    annotate("text", x = 5.1e-7, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(3)", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range_slab, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_gradient2(low = "#A50026", mid = "white", high = "#313695", limits = displacement_range_slab) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = "Displacement (km)") +
    theme_bw(base_size = 14) +
    theme_2() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank())

  p2 <- df_slab_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = width)) +
    geom_rect() +
    rects_slab_b +
    annotate("text", x = 1.1e-1, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(1)", fontface = "bold") +
    annotate("text", x = 1.3e-4, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(2)", fontface = "bold") +
    annotate("text", x = 5.1e-7, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(3)", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range_slab, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(option = "mako", limits = width_range_slab, direction = -1) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = "Width (km)") +
    theme_bw(base_size = 14) +
    theme_2()

  p3 <- df_slab_tiled |>
    ggplot(aes(xmin = xmin, xmax = xmax, ymin = B_factor - 1.0, ymax = B_factor + 1.0, fill = max_velocity)) +
    geom_rect() +
    rects_slab_w +
    annotate("text", x = 1.1e-1, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(1)", color = "white", fontface = "bold") +
    annotate("text", x = 1.3e-4, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(2)", color = "white", fontface = "bold") +
    annotate("text", x = 5.1e-7, y = 4, size = 4, hjust = 0.5, vjust = 0.5, label = "(3)", color = "white", fontface = "bold") +
    scale_x_continuous(trans = "log10", labels = label_log(), breaks = breaks_z, limits = z_range_slab, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
    annotation_logticks(sides = "b", linewidth = 0.2) +
    scale_fill_viridis_c(option = "rocket", limits = max_velocity_range_slab) +
    labs(x = bquote("Z (" * s^-1 * K^-1 * ")"), y = "B", fill = bquote("Max " * italic(u)[y] * " (cm/yr)")) +
    theme_bw(base_size = 14) +
    theme_2() +
    theme(axis.text = element_blank(), axis.title = element_blank(), legend.title = element_text(margin = margin(0, 0, 2, 0)))

  out_path_comp <- str_replace(out_path, ".png", "-comp.png")

  if (plot_exists(out_path_comp)) {
    return(invisible())
  }

  p <- (p0 | p3) / (p2 | p1)
  suppressWarnings(ggsave(out_path_comp, plot = p, width = 4.5, height = 5.0, dpi = 300, bg = "white"))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_1 <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey90"),
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.tag.location = "panel",
    plot.tag.position = "topleft",
    plot.tag = element_text(size = 18, margin = margin(5, 0, 0, 0), hjust = 0, color = "black", face = "bold"),
    axis.ticks = element_blank(),
    legend.justification = "left",
    legend.position = "inside",
    legend.position.inside = c(0.07, 0.75),
    legend.direction = "horizontal",
    legend.key.height = unit(0.6, "lines"),
    legend.key.width = unit(1.2, "lines"),
    legend.ticks = element_line(color = "black", linewidth = 0.4),
    legend.ticks.length = unit(0.1, "lines"),
    legend.frame = element_rect(color = "black", linewidth = 0.4),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.box.spacing = unit(0.1, "lines"),
    legend.box.margin = margin(),
    legend.margin = margin(-13, 0, 0, 0),
    legend.title = element_text(hjust = 0, vjust = 0, size = 12, margin = margin(0, 0, 5, 0)),
    legend.title.position = "top",
    legend.text = element_text(size = 11, margin = margin(2, 0, 0, 0)),
    legend.background = element_blank()
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_2 <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey90"),
    plot.margin = margin(2, 5, 2, 5),
    plot.title = element_text(siz = 14, hjust = 0.5),
    plot.tag.location = "panel",
    plot.tag.position = "topleft",
    plot.tag = element_text(size = 18, margin = margin(5, 0, 0, 0), hjust = 0, color = "black", face = "bold"),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.height = unit(0.6, "lines"),
    legend.key.width = unit(1.2, "lines"),
    legend.ticks = element_line(color = "black", linewidth = 0.4),
    legend.ticks.length = unit(0.1, "lines"),
    legend.frame = element_rect(color = "black", linewidth = 0.4),
    legend.box.spacing = unit(0.1, "lines"),
    legend.box.margin = margin(),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(hjust = 0.5, vjust = 0, size = 10, margin = margin(0, 0, 5, 0)),
    legend.title.position = "top",
    legend.text = element_text(size = 11, margin = margin(2, 0, 0, 0)),
    legend.background = element_blank()
  )
}
