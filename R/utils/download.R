#######################################################
## Download ASPECT results from OSF repo         !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
download_simulation_results_from_osf <- function(out_dir_csv, out_dir_sim) {
  if (!file.exists(file.path(out_dir_csv, "depth-profile-data.csv"))) {
    if (!dir.exists(out_dir_csv)) dir.create(out_dir_csv, recursive = TRUE)
    osf_retrieve_node("9phwc") |>
      osf_ls_files(pattern = "csv") |>
      osf_download(path = out_dir_csv, conflicts = "overwrite")
  } else {
    cat(" -- Found data: ", file.path(out_dir_csv, "depth-profile-data.csv"), "\n", sep = "")
  }

  if (!dir.exists(out_dir_sim)) dir.create(out_dir_sim, recursive = TRUE)

  patterns <- c(
    "Z3.0e0_B2", "Z4.7e2_B2", "Z7.0e7_B2",
    "Z3.0e0_B4", "Z4.7e2_B4", "Z7.0e7_B4",
    "Z3.0e0_B6", "Z4.7e2_B6", "Z7.0e7_B6",
    "Z3.0e0_B8", "Z4.7e2_B8", "Z7.0e7_B8",
    "Z3.0e0_B10", "Z4.7e2_B10", "Z7.0e7_B10"
  )

  walk(patterns, ~ {
    plume_dir <- file.path(out_dir_sim, paste0("plume_", .x))
    slab_dir <- file.path(out_dir_sim, paste0("slab_", .x))

    need_download <- !(dir.exists(plume_dir) && dir.exists(slab_dir))

    if (need_download) {
      osf_retrieve_node("9phwc") |>
        osf_ls_files(pattern = "results") |>
        osf_ls_files(pattern = .x) |>
        osf_download(path = out_dir_sim, conflicts = "overwrite", recurse = TRUE)
    } else {
      cat(" -- Found data: ", plume_dir, "\n", sep = "")
      cat(" -- Found data: ", slab_dir, "\n", sep = "")
    }
  })
}
