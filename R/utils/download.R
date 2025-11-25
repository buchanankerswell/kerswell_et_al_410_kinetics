#######################################################
## Download ASPECT results from OSF repo         !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
download_simulation_results_from_osf <- function(out_dir) {
  if (!file.exists(paste0(out_dir, "/depth-profile-data.csv"))) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    osf_retrieve_node("9phwc") |>
      osf_ls_files(pattern = "csv") |>
      osf_download(path = out_dir, conflicts = "overwrite")
  } else {
    cat(" -- Found data: ", paste0(out_dir, "/depth-profile-data.csv"), "\n", sep = "")
  }
}
