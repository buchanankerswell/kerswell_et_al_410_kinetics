#!/usr/bin/env Rscript

#######################################################
## Visualize 410 structure                       !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 3) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript simulation.R [util_dir] [data_dir] [out_dir]\n")
    return(invisible(NULL))
  }

  util_dir <- args[1]
  data_dir <- args[2]
  out_dir <- args[3]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  in_data <- file.path(data_dir, "depth-profile-data.csv")
  out_path <- file.path(out_dir, "410-structure.png")

  cat("    --------------------------------------------------\n")
  cat("    Drawing 410 structure summary\n")
  cat("    --------------------------------------------------\n")

  visualize_410_structure(in_data, out_path)
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
