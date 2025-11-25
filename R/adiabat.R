#!/usr/bin/env Rscript

#######################################################
## Visualize Adiabatic Reference Conditions      !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 4) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript adiabat.R [util_dir] [model_id] [data_dir] [out_dir]\n")
    return(invisible(NULL))
  }

  util_dir <- args[1]
  model_id <- args[2]
  data_dir <- args[3]
  out_dir <- args[4]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  in_adiabatic_profile <- file.path(data_dir, paste0(model_id, "-adiabatic-profile.txt"))
  in_material_table <- file.path(data_dir, paste0(model_id, "-material-table.txt"))
  in_reaction_data <- file.path(data_dir, "reaction-410-profile.txt")

  out_material_table <- file.path(out_dir, paste0(model_id, "-material-table.png"))
  out_material_profile <- file.path(out_dir, "material-property-profile.png")
  out_thermodynamic_profile <- file.path(out_dir, "thermodynamic-property-profile.png")

  missing <- c()
  if (!file.exists(in_material_table)) missing <- c(missing, in_material_table)
  if (!file.exists(in_reaction_data)) missing <- c(missing, in_reaction_data)

  if (length(missing) > 0) {
    cat("    --------------------------------------------------\n")
    cat(" !! Warning: the following input files do not exist:\n")
    for (f in missing) cat(" -- ", f, "\n", sep = "")
    return(invisible(NULL))
  }

  cat("    --------------------------------------------------\n")
  cat("    Drawing reference material properties\n")
  cat("    --------------------------------------------------\n")

  tryCatch(
    {
      visualize_material_table(in_adiabatic_profile, in_material_table, out_material_table)
      visualize_material_profile(in_reaction_data, out_material_profile)
      visualize_thermodynamic_profile(in_reaction_data, out_thermodynamic_profile)
    },
    error = function(e) {
      cat("    --------------------------------------------------\n")
      cat(" !! Error: drawing issue: ", conditionMessage(e), "\n", sep = "")
    }
  )
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
