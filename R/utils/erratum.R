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

  # Water content (ppm) on a contour of constant Z at fixed Gamma
  coh_on_contour <- function(d_m, z, gamma) {
    (z * d_m / (6.67 * gamma))^(1 / n)
  }

  # Z ranges: published (simulation/Makefile Z_FACTOR endpoints) and
  # effective (scaled). Permitted edges for the extended dry floor:
  # low edge from 10-mm grains and 1 ppm OH at -2 sigma,
  # high edge from 1-mm grains and 5000 ppm at +2 sigma.
  z_pub <- c(3.0e0, 7.0e7)
  z_eff <- z_pub / sec_per_year
  z_env <- c(
    (6.67 / 1e-2) * exp(ln_gamma_low) * 1^n,
    (6.67 / 1e-3) * exp(ln_gamma_high) * 5000^n
  )

  # Facet Gammas: -2 sigma, best fit, +2 sigma (self-consistent per panel)
  facet_spec <- tribble(
    ~ln_gamma, ~facet,
    ln_gamma_low, "ln~Gamma==-25.6~group('(', -2*sigma, ')')",
    ln_gamma_best, "ln~Gamma==-18.0~(best~fit)",
    ln_gamma_high, "ln~Gamma==-10.4~group('(', +2*sigma, ')')"
  ) |>
    mutate(facet = factor(facet, levels = facet))

  # Log10 shift of Gamma relative to best fit: contours move with Gamma so
  # each panel keeps the same count of on-scale power-of-10 contours
  gamma_shift <- function(ln_gamma) round((ln_gamma - ln_gamma_best) / log(10))

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

  band_levels <- c("published", "effective")

  df_bands <- expand_grid(facet_spec, band = band_levels) |>
    mutate(
      band = factor(band, levels = band_levels),
      facet = factor(facet, levels = facet_spec$facet),
      z_lo = if_else(band == "published", z_pub[1], z_eff[1]),
      z_hi = if_else(band == "published", z_pub[2], z_eff[2])
    ) |>
    pmap_dfr(function(ln_gamma, facet, band, z_lo, z_hi, ...) {
      band_quad(ln_gamma, z_lo, z_hi) |>
        mutate(facet = facet, band = band)
    }) |>
    mutate(facet = factor(facet, levels = facet_spec$facet))

  # Self-consistent Z contours per facet: powers of 10 shifted by the
  # Gamma offset so each panel shows the same number of on-scale contours
  df_contours <- facet_spec |>
    mutate(shift = gamma_shift(ln_gamma)) |>
    pmap_dfr(function(ln_gamma, facet, shift, ...) {
      expand_grid(
        z = 10^seq(-9 + shift, 12 + shift, by = 3),
        d_m = d_range
      ) |>
        mutate(
          d = d_m * 1e3, # display grain size in mm
          COH = coh_on_contour(d_m, z, exp(ln_gamma)),
          facet = facet
        )
    }) |>
    mutate(facet = factor(facet, levels = facet_spec$facet))

  # Contour labels just inside the right edge of each panel
  df_contour_labels <- facet_spec |>
    mutate(shift = gamma_shift(ln_gamma)) |>
    pmap_dfr(function(ln_gamma, facet, shift, ...) {
      tibble(z = 10^seq(-9 + shift, 12 + shift, by = 3)) |>
        mutate(
          d = (d_range[2] * 1e3) / 10^0.15,
          d_m = d / 1e3,
          COH = coh_on_contour(d_m, z, exp(ln_gamma)),
          label = paste0("10^", log10(z)),
          facet = facet
        )
    }) |>
    mutate(facet = factor(facet, levels = facet_spec$facet)) |>
    filter(COH > coh_range[1] * 1.5, COH < coh_range[2] / 1.5)

  # Permitted limit contours, one per originating facet: low edge at -2 sigma,
  # high edge at +2 sigma. Each spans the full d range at its own Gamma.
  df_limits <- tribble(
    ~ln_gamma, ~facet_id, ~z, ~label,
    ln_gamma_low, 1L, z_env[1], "5.1%*%10^-9",
    ln_gamma_high, 3L, z_env[2], "1.4%*%10^11"
  ) |>
    mutate(
      facet = factor(facet_spec$facet[facet_id], levels = facet_spec$facet)
    ) |>
    pmap_dfr(function(ln_gamma, facet, z, label, ...) {
      tibble(d_m = seq(log10(d_range[1]), log10(d_range[2]), length.out = 100)) |>
        mutate(
          d_m = 10^d_m,
          d = d_m * 1e3,
          COH = coh_on_contour(d_m, z, exp(ln_gamma)),
          z = z,
          label = label,
          facet = facet
        )
    }) |>
    mutate(facet = factor(facet, levels = facet_spec$facet)) |>
    mutate(boundary = "permitted")

  # Drop background labels that collide with the limit labels:
  # 1e-9 beside 5.1e-09 (left), 1e12 beside 1.4e11 (right)
  df_contour_labels <- df_contour_labels |>
    filter(
      !(facet == facet_spec$facet[1] & z == 1e-9),
      !(facet == facet_spec$facet[3] & z == 1e12)
    )

  # Labels right of center, same baseline style as the background contour labels
  # but shifted left so the longer limit strings clear the panel edge
  df_limit_labels <- tribble(
    ~ln_gamma, ~facet_id, ~z, ~label,
    ln_gamma_low, 1L, z_env[1], "5.1%*%10^-9",
    ln_gamma_high, 3L, z_env[2], "1.4%*%10^11"
  ) |>
    mutate(
      facet = factor(facet_spec$facet[facet_id], levels = facet_spec$facet),
      d = (d_range[2] * 1e3) / 10^0.55,
      d_m = d / 1e3,
      COH = coh_on_contour(d_m, z, exp(ln_gamma))
    )

  # Condition boxes (independent of Gamma), drawn identically in every facet:
  # actual Hosoya/manuscript box (d = 1-10 mm, COH = 50-5000 ppm) and
  # extended dry box (d = 1-10 mm, COH = 1-5000 ppm)
  df_experiment <- tribble(
    ~d, ~COH, ~boundary,
    1, 50, "hosoya", 10, 50, "hosoya", 10, 5000, "hosoya",
    1, 5000, "hosoya", 1, 50, "hosoya",
    1, 1, "extended", 10, 1, "extended", 10, 5000, "extended",
    1, 5000, "extended", 1, 1, "extended"
  )

  # Downshift arrow in the best-fit facet only: the published Z = 3.0e0
  # corner at d = 10 mm ran as Z' = 9.5e-8, a factor 3.15e7 lower in Z
  # (220x lower in COH)
  df_arrow <- tibble(
    x = 0.8,
    y_top = coh_on_contour(1e-2, z_pub[1], exp(ln_gamma_best)),
    y_bot = coh_on_contour(1e-2, z_eff[1], exp(ln_gamma_best)),
    facet = factor(facet_spec$facet[2], levels = facet_spec$facet)
  ) |>
    mutate(y_mid = sqrt(y_top * y_bot))

  # Arrow note in the best-fit facet
  df_notes <- tibble(
    x = 0.7,
    y = df_arrow$y_mid,
    vjust = 0.5,
    hjust = 1,
    label = "Accidental\nscaling\n(red to blue)",
    facet = df_arrow$facet
  )

  p <- ggplot(df_bands, aes(x = d, y = COH)) +
    geom_line(
      data = df_contours,
      aes(x = d, y = COH, group = z),
      color = "grey70",
      linewidth = 0.3,
      linetype = "dashed"
    ) +
    geom_text(
      data = df_contour_labels,
      aes(x = d, y = COH, label = label),
      parse = TRUE,
      hjust = 0.5,
      vjust = 0,
      size = 3,
      color = "grey30"
    ) +
    geom_polygon(
      data = df_bands,
      aes(group = band, fill = band),
      color = NA,
      alpha = 0.5,
      linewidth = 0.4
    ) +
    geom_path(
      data = filter(df_experiment, boundary == "extended"),
      aes(x = d, y = COH, linetype = boundary, group = boundary),
      color = "black",
      linewidth = 0.7
    ) +
    geom_path(
      data = filter(df_experiment, boundary == "hosoya"),
      aes(x = d, y = COH, linetype = boundary, group = boundary),
      color = "black",
      linewidth = 0.5
    ) +
    geom_line(
      data = df_limits,
      aes(x = d, y = COH, linetype = boundary, group = facet),
      color = "black",
      linewidth = 0.9
    ) +
    geom_text(
      data = df_limit_labels,
      aes(x = d, y = COH, label = label),
      parse = TRUE,
      hjust = 0.5,
      vjust = -0.5,
      size = 3,
      color = "black",
      fontface = "bold"
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
    facet_wrap(vars(facet), nrow = 1, labeller = label_parsed) +
    scale_x_log10(
      breaks = 10^(-0:1),
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
      name = NULL,
      breaks = band_levels,
      labels = list(
        expression("Published " * Z),
        expression("Effective " * Z * minute)
      ),
      values = setNames(
        c("#E41A1C", "#377EB8"),
        band_levels
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      breaks = c("permitted", "hosoya", "extended"),
      labels = list(
        expression("Permitted " * Z * " range (" * "" %+-% 2 * sigma * ")"),
        "Hosoya experiments",
        "Extended dry conditions"
      ),
      values = c(
        "permitted" = "solid", "hosoya" = "solid", "extended" = "dashed"
      )
    ) +
    labs(
      x = bquote("Grain size " * italic(d) * " (mm)"),
      y = bquote(C["OH"] * " (ppm)")
    ) +
    guides(
      fill = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(
          color = "black",
          linewidth = c(0.9, 0.5, 0.7)
        )
      )
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.width = grid::unit(1.5, "lines"),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.box.margin = margin(),
      legend.margin = margin(),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text = element_text(size = 13)
    )

  ggsave(out_path, plot = p, width = 11, height = 4.2, dpi = 300, bg = "white")
}
