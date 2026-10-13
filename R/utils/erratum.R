#######################################################
## Visualize Erratum Parameter Ranges            !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_erratum_parameter_ranges <- function(out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  # Kinetic parameters (manuscript Section 2.2.3.2 and erratum)
  n <- 3.2
  sec_per_year <- 3.15e7 # s/yr scaling behind the effective Z' values
  ln_gamma_best <- -18.0 # best-fit growth-rate prefactor
  ln_gamma_low <- -25.6 # -2 sigma
  ln_gamma_high <- -10.4 # +2 sigma

  # Axis ranges
  d_range <- c(1e-4, 1e-1) # grain size (m), i.e. 0.1 to 100 mm
  coh_range <- c(1e-2, 10^5.5) # water content (ppm)

  # Z = (6.67 / d) * Gamma * COH^n
  kinetic_prefactor <- function(d_m, coh, ln_gamma) {
    (6.67 / d_m) * exp(ln_gamma) * coh^n
  }

  # Water content (ppm) on a contour of constant Z at fixed Gamma
  coh_on_contour <- function(d_m, z, gamma) {
    (z * d_m / (6.67 * gamma))^(1 / n)
  }

  # Z ranges: published (simulation/Makefile Z_FACTOR endpoints),
  # effective (scaled), and full permitted envelope
  z_pub <- c(3.0e0, 7.0e7)
  z_eff <- z_pub / sec_per_year
  z_env <- c(
    kinetic_prefactor(1e-2, 1, ln_gamma_low), # 10-mm grains, 1 ppm OH
    kinetic_prefactor(1e-3, 5000, ln_gamma_high) # 1-mm grains, 5000 ppm OH
  )

  band_spec <- tribble(
    ~band, ~z_lo, ~z_hi, ~ln_gamma, ~style,
    "envelope", z_env[1], z_env[2], ln_gamma_best, "envelope",
    "published", z_pub[1], z_pub[2], ln_gamma_best, "solid",
    "effective", z_eff[1], z_eff[2], ln_gamma_best, "solid",
    "published_ghost", z_pub[1], z_pub[2], ln_gamma_high, "ghost",
    "effective_ghost", z_eff[1], z_eff[2], ln_gamma_low, "ghost"
  )

  band_levels <- band_spec$band # legend/draw order

  # (d, COH) quad for Z in [z_lo, z_hi] at fixed Gamma over d in [d_lo, d_hi]:
  # bottom-left, bottom-right, top-right, top-left corners in order
  band_quad <- function(ln_gamma, z_lo, z_hi, d_lo = 1e-3, d_hi = 1e-2) {
    gamma <- exp(ln_gamma)
    tibble(
      d_m = c(d_lo, d_hi, d_hi, d_lo),
      z = c(z_lo, z_lo, z_hi, z_hi)
    ) |>
      mutate(
        d = d_m * 1e3, # display grain size in mm
        COH = coh_on_contour(d_m, z, gamma)
      )
  }

  df_bands <- band_spec |>
    mutate(band = factor(band, levels = band_levels)) |>
    pmap_dfr(function(band, z_lo, z_hi, ln_gamma, style, ...) {
      band_quad(ln_gamma, z_lo, z_hi) |>
        mutate(band = band, style = style)
    })

  # Best-fit Z contours every 10^3 (K^-1 s^-1), drawn first as background texture
  df_best <- expand_grid(
    z = 10^seq(-9, 12, by = 3),
    d_m = d_range
  ) |>
    mutate(
      d = d_m * 1e3, # display grain size in mm
      COH = coh_on_contour(d_m, z, exp(ln_gamma_best))
    )

  # Contour labels just inside the right edge of the panel
  df_best_labels <- tibble(z = 10^seq(-9, 12, by = 3)) |>
    mutate(
      d = (d_range[2] * 1e3) / 10^0.15,
      d_m = d / 1e3,
      COH = coh_on_contour(d_m, z, exp(ln_gamma_best)),
      label = paste0("10^", log10(z))
    ) |>
    filter(COH > coh_range[1] * 1.5, COH < coh_range[2] / 1.5)

  # Downshift arrow: the published Z = 3.0e0 corner at d = 10 mm ran as
  # Z' = 9.5e-8, a factor 3.15e7 lower in Z (221x lower in COH)
  df_arrow <- tibble(
    x = 0.8,
    y_top = coh_on_contour(1e-2, z_pub[1], exp(ln_gamma_best)),
    y_bot = coh_on_contour(1e-2, z_eff[1], exp(ln_gamma_best))
  ) |>
    mutate(y_mid = sqrt(y_top * y_bot))

  # Margin notes in the empty ground on either side of the bands
  df_notes <- tribble(
    ~x, ~y, ~vjust, ~hjust, ~label,
    0.7, df_arrow$y_mid, 0.5, 1, "Accidental\nscaling",
    0.14, 90000, 1, 0, "Permitted range\nHosoya (2005)"
  )

  # Overlap note at the center of the overlapping ghost shapes
  df_overlap <- tibble(
    x = 1.1e1,
    y = 50,
    label = "Bands\noverlap\nat 95% CI"
  )

  p <- ggplot(df_bands, aes(x = d, y = COH)) +
    geom_line(
      data = df_best,
      aes(x = d, y = COH, group = z),
      color = "grey80",
      linewidth = 0.3,
      linetype = "dashed"
    ) +
    geom_text(
      data = df_best_labels,
      aes(x = d, y = COH, label = label),
      parse = TRUE,
      hjust = 0.5,
      vjust = 0,
      size = 3,
      color = "grey30"
    ) +
    geom_polygon(
      data = filter(df_bands, style == "envelope"),
      aes(group = band, fill = band),
      color = "black",
      alpha = 0.50,
      linewidth = 0.4
    ) +
    geom_polygon(
      data = filter(df_bands, style == "ghost"),
      aes(group = band, fill = band),
      color = NA,
      alpha = 0
    ) +
    geom_path(
      data = filter(df_bands, band == "published_ghost"),
      aes(group = band),
      color = "#E41A1C",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_path(
      data = filter(df_bands, band == "effective_ghost"),
      aes(group = band),
      color = "#377EB8",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_polygon(
      data = filter(df_bands, style == "solid"),
      aes(group = band, fill = band),
      color = "black",
      alpha = 0.5,
      linewidth = 0.4
    ) +
    geom_segment(
      data = df_arrow,
      aes(x = x, xend = x, y = y_top, yend = y_bot),
      color = "grey30",
      linewidth = 0.4,
      arrow = grid::arrow(length = grid::unit(0.25, "cm"), type = "closed")
    ) +
    geom_text(
      data = df_notes,
      aes(x = x, y = y, label = label, vjust = vjust, hjust = hjust),
      size = 3,
      color = "grey30"
    ) +
    geom_text(
      data = df_overlap,
      aes(x = x, y = y, label = label),
      hjust = 0,
      vjust = 0.5,
      size = 3,
      color = "black",
    ) +
    scale_x_log10(
      breaks = 10^(-1:2),
      labels = label_log()
    ) +
    scale_y_log10(
      breaks = 10^(-2:5),
      labels = label_log()
    ) +
    coord_cartesian(
      xlim = d_range * 1e3,
      ylim = coh_range,
      expand = FALSE
    ) +
    annotation_logticks(sides = "bl", linewidth = 0.2) +
    scale_fill_manual(
      name = "Z range",
      breaks = band_levels,
      labels = list(
        expression("Uncertainty envelope (" %+-% 2 * sigma * ")"),
        expression("Published " * Z),
        expression("Effective " * Z * minute),
        expression("Published " * Z * " (" * +2 * sigma * ")"),
        expression("Effective " * Z * minute * " (" * -2 * sigma * ")")
      ),
      values = setNames(
        c("grey60", "#E41A1C", "#377EB8", "#E41A1C", "#377EB8"),
        band_levels
      )
    ) +
    labs(
      x = bquote("Grain size " * italic(d) * " (mm)"),
      y = bquote(C["OH"] * " (ppm)")
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.background = element_blank(),
      legend.title.position = "top"
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          colour = c("black", "black", "black", "#E41A1C", "#377EB8"),
          linetype = c("solid", "solid", "solid", "dashed", "dashed")
        )
      )
    )

  ggsave(out_path, plot = p, width = 6.5, height = 3.5, dpi = 300, bg = "white")
}
