#######################################################
## General Helper Functions                      !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ensure_output_dir <- function(out_path) {
  parent_dir <- dirname(out_path)
  if (!dir.exists(parent_dir)) {
    dir.create(parent_dir)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_exists <- function(out_path) {
  if (file.exists(out_path)) {
    cat(" -- Found plot: ", basename(out_path), "!\n", sep = "")
    return(TRUE)
  }

  ensure_output_dir(out_path)

  cat(" -> ", basename(out_path), "\n", sep = "")
  FALSE
}
