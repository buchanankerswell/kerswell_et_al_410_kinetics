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
source(file.path(this_dir, "statistics.R"))
source(file.path(this_dir, "viscosity.R"))
source(file.path(this_dir, "averages.R"))

#######################################################
## .1. Visualize All Simulation Plots            !!! ##
#######################################################
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
  prefix <-
    str_replace_all(basename(in_dir), "_", "-") |>
    str_replace("output-", "")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  in_stat <- file.path(in_dir, "statistics")
  in_res <- file.path(in_dir, "stokes_residuals.txt")
  in_davg <- file.path(in_dir, "depth_average.txt")

  out_stats <- file.path(out_dir, paste0(prefix, "-statistics.png"))
  out_viscosity <- file.path(out_dir, paste0(prefix, "-viscosity-profile.png"))
  out_depth <- file.path(out_dir, paste0(prefix, "-depth-averages.png"))

  missing <- c()
  # if (!file.exists(in_stat)) missing <- c(missing, in_stat)
  # if (!file.exists(in_res)) missing <- c(missing, in_res)
  if (!file.exists(in_davg)) missing <- c(missing, in_davg)

  if (length(missing) > 0) {
    cat("    --------------------------------------------------\n")
    cat(" !! Warning: the following input files do not exist:\n")
    for (f in missing) {
      cat(" -- ", f, "\n", sep = "")
    }
    return(invisible(NULL))
  }

  cat("    --------------------------------------------------\n")
  cat("    Processing simulation: ", basename(in_dir), "\n", sep = "")
  cat("    --------------------------------------------------\n")

  tryCatch(
    {
      # visualize_statistics(in_stat, in_res, out_stats)
      visualize_viscosity_profile(in_davg, out_viscosity)
      visualize_depth_averages(in_davg, out_depth)
    },
    error = function(e) {
      cat("    --------------------------------------------------\n")
      cat(" !! Error: drawing issue: ", conditionMessage(e), "\n", sep = "")
    }
  )
}

if (
  !interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))
) {
  main()
}
