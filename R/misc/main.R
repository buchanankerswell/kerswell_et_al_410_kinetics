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
source(file.path(this_dir, "deformation.R"))
source(file.path(this_dir, "combined-viscosity.R"))

#######################################################
## .1. Visualize All Perplex/Burnman Plots       !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript main.R [in_dir] [out_dir]\n")
    cat(" !! Example: Rscript main.R /path/to/sim_out_dir /path/to/fig_dir\n")
    return(invisible(NULL))
  }

  profile_paths <- args[1:(length(args) - 2)]
  out_dir_profiles <- args[length(args) - 1]
  out_dir_deformation <- args[length(args)]

  if (!dir.exists(out_dir_profiles)) {
    dir.create(out_dir_profiles, recursive = TRUE)
  }

  if (!dir.exists(out_dir_deformation)) {
    dir.create(out_dir_deformation, recursive = TRUE)
  }

  out_deformation <- file.path(out_dir_deformation, "deformation-map.png")
  out_profiles <- file.path(out_dir_profiles, "combined-viscosity-profiles.png")

  cat("    --------------------------------------------------\n")
  cat("    Processing misc plots\n")
  cat("    --------------------------------------------------\n")

  tryCatch(
    {
      visualize_deformation_map(out_deformation)
      visualize_viscosity_profiles(profile_paths, out_profiles)
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
