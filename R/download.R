#!/usr/bin/env Rscript

#######################################################
## Download ASPECT results from OSF repo         !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript main.R [util_dir] [out_dir]\n")
    return(invisible(NULL))
  }

  util_dir <- args[1]
  out_dir <- args[2]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  cat("    --------------------------------------------------\n")
  cat("    Downloading simulation results from OSF repo\n")
  cat("    --------------------------------------------------\n")

  download_simulation_results_from_osf(out_dir)
}

if (
  !interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))
) {
  main()
}
