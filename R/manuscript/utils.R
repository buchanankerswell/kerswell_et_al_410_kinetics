#!/usr/bin/env Rscript

#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(knitr)
})

#######################################################
## .1. Utility Functions                         !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
convert_preexponential <- function(df) {
  df |>
    mutate(A = case_when(
      !is.na(A_mpa_coef) & !is.na(A_mpa_exp) ~ {
        A_mpa_coef * 10^A_mpa_exp * 10^(-(6 * n)) * 10^(-(6 * m))
      },
      !is.na(A_pa_coef) & !is.na(A_pa_exp) ~ {
        A_pa_coef * 10^A_pa_exp
      },
      TRUE ~ NA_real_
    ), .after = type) |>
    select(-all_of(c("A_mpa_coef", "A_mpa_exp", "A_pa_coef", "A_pa_exp")))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_markdown_table <- function(df, out_file) {
  df |>
    mutate(across(everything(), ~ ifelse(is.na(.), "-", .))) |>
    mutate(A = format(A, scientific = TRUE, digits = 3)) |>
    mutate(Ea = format(Ea, scientific = TRUE, digits = 3)) |>
    kable(format = "markdown") |>
    writeLines(out_file)
}
