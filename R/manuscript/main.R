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
## .1. Write Rheological Parameters to md Table  !!! ##
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

  in_path <- args[1]
  out_path <- args[2]

  if (!file.exists(in_path)) {
    cat("    --------------------------------------------------\n")
    cat(" !! Warning: the following input files do not exist:\n")
    cat(" -- ", in_path, "\n", sep = "")
    return(invisible(NULL))
  }

  cat("    --------------------------------------------------\n")
  cat("    Processing csv file: ", basename(in_path), "\n", sep = "")

  tryCatch(
    {
      read_csv(in_path, na = "NA", show_col_type = FALSE) |>
        convert_preexponential() |>
        write_markdown_table(out_path)
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
