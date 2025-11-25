#######################################################
## Download ASPECT results from OSF repo         !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
download_simulation_results_from_osf <- function(out_dir) {
  osf_retrieve_node("9phwc") |>
    osf_ls_files(project, pattern = "results") |>
    osf_download(path = out_dir, recurse = TRUE, conflicts = "overwrite", progress = TRUE)

}
