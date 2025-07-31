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
## .1. Visualize Deflection                      !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visualize_deflection <- function(deflection_path, out_path) {
  if (plot_exists(out_path)) {
    return(invisible())
  }

  df <-
    read_deflections(deflection_path) |>
    filter(timestep == 100) |>
    mutate(model_id = str_replace_all(model_id, "_", "-"))
  df_slab <- df |> filter(str_detect(model_id, "slab"))
  df_plume <- df |> filter(str_detect(model_id, "plume"))

  p0 <-
    df_slab |>
    ggplot(aes(x = deflection, y = sharpness, fill = model_id)) +
    geom_point(shape = 21, size = 2.5, stroke = 0.5, color = "black") +
    labs(
      x = "Distance from Equil. Depth (km)",
      y = "Phase Transition Zone Width (km)",
      fill = NULL,
      title = "Slabs"
    ) +
    theme_bw(base_size = 14) +
    theme_2()

  p1 <-
    df_plume |>
    ggplot(aes(x = deflection, y = sharpness, fill = model_id)) +
    geom_point(shape = 21, size = 2.5, stroke = 0.5, color = "black") +
    labs(
      x = "Distance from Equil. Depth (km)",
      y = "Phase Transition Zone Width (km)",
      fill = NULL,
      title = "Plumes"
    ) +
    theme_bw(base_size = 14) +
    theme_2()

  p <- p0 / p1

  ggsave(
    out_path,
    plot = p,
    width = 8.5,
    height = 8.5,
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

  in_deflection <- file.path(in_dir, "deflection_sharpness_data.csv")
  out_deflection <- file.path(out_dir, "deflection-sharpness.png")
  visualize_deflection(in_deflection, out_deflection)
}

if (
  !interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))
) {
  main()
}
