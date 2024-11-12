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
    stop("Cannot determine script location.")
  }
}

this_dir <- get_script_dir()

source(file.path(this_dir, "utils.R"))
source(file.path(this_dir, "profile.R"))
source(file.path(this_dir, "table.R"))

#######################################################
## .1. Visualize All Perplex/Burnman Plots       !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    cat("  Usage: Rscript main.R [in_dir] [out_dir]\n")
    cat("  Example: Rscript main.R /path/to/sim_out_dir /path/to/fig_dir\n")
    return(invisible(NULL))
  }

  model_id <- args[1]
  in_dir <- args[2]
  out_dir <- args[3]

  if (!dir.exists(out_dir)) {
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    dir.create(out_dir, recursive = TRUE)
    cat("==> Created out directory:", out_dir, "\n")
  }

  in_adiabatic_profile <-
    file.path(in_dir, paste0(model_id, "-adiabatic-profile.txt"))
  in_material_table <-
    file.path(in_dir, paste0(model_id, "-material-table.txt"))
  in_driving_force_profile <-
    file.path(in_dir, "thermodynamic-driving-force-profile-olivine-wadsleyite.txt")

  out_adiabatic_profile <-
    file.path(out_dir, paste0(model_id, "-adiabatic-profile.png"))
  out_material_table <-
    file.path(out_dir, paste0(model_id, "-material-table.png"))
  out_driving_force_profile <-
    file.path(out_dir, paste0(model_id, "-driving-force-profile.png"))

  missing <- c()
  if (!file.exists(in_adiabatic_profile)) {
    missing <- c(missing, in_adiabatic_profile)
  }
  if (!file.exists(in_material_table)) {
    missing <- c(missing, in_material_table)
  }
  if (!file.exists(in_driving_force_profile)) {
    missing <- c(missing, in_driving_force_profile)
  }

  if (length(missing) > 0) {
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    cat("Warning: The following input files do not exist:\n")
    for (f in missing) {
      cat(" - ", f, "\n")
    }
    cat("Please check the paths and try again.\n")
    return(invisible(NULL))
  }

  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("==> Processing perplex model:", basename(in_dir), "\n")

  tryCatch(
    {
      visualize_adiabatic_profile(
        in_adiabatic_profile,
        out_adiabatic_profile
      )
      visualize_material_table(
        in_adiabatic_profile,
        in_material_table,
        out_material_table
      )
      visualize_driving_force_profile(
        in_driving_force_profile,
        out_driving_force_profile
      )
    },
    error = function(e) {
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
      cat("  Error occurred during visualization:", conditionMessage(e), "\n")
    }
  )
}

if (
  !interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))
) {
  main()
}
