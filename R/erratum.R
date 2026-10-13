#!/usr/bin/env Rscript

#######################################################
## Visualize Erratum Parameter Ranges            !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript erratum.R [util_dir] [out_dir]\n")
    return(invisible(NULL))
  }

  util_dir <- args[1]
  out_dir <- args[2]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_path <- file.path(out_dir, "erratum-kinetics.png")

  cat("    --------------------------------------------------\n")
  cat("    Drawing erratum parameter ranges\n")
  cat("    --------------------------------------------------\n")

  visualize_erratum_parameter_ranges(out_path)
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
