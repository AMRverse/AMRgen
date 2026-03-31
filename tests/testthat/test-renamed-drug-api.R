# ===================================================================== #
#  Licensed as GPL-v3.0.                                                #
#                                                                       #
#  Developed as part of the AMRverse (https://github.com/AMRverse):     #
#  https://github.com/AMRverse/AMRgen                                   #
#                                                                       #
#  We created this package for both routine data analysis and academic  #
#  research and it was publicly released in the hope that it will be    #
#  useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                       #
#  This R package is free software; you can freely use and distribute   #
#  it for both personal and commercial purposes under the terms of the  #
#  GNU General Public License version 3.0 (GNU GPL-3), as published by  #
#  the Free Software Foundation.                                        #
# ===================================================================== #

test_that("get_binary_matrix supports drug, pheno_drug, and geno_class", {
  pheno_table <- tibble::tibble(
    id = c("s1", "s2"),
    drug = AMR::as.ab(c("Ciprofloxacin", "Ciprofloxacin")),
    pheno_clsi = AMR::as.sir(c("R", "S")),
    ecoff = AMR::as.sir(c("NWT", "WT")),
    mic = AMR::as.mic(c("4", "0.25"))
  )

  geno_table <- tibble::tibble(
    id = c("s1", "s2"),
    marker = c("gyrA_S83F", "parC_S80I"),
    drug_class = c("Quinolones", "Quinolones"),
    drug = AMR::as.ab(c("Ciprofloxacin", "Ciprofloxacin"))
  )

  result <- get_binary_matrix(
    geno_table = geno_table,
    pheno_table = pheno_table,
    pheno_drug = "Ciprofloxacin",
    geno_class = "Quinolones",
    sir_col = "pheno_clsi",
    keep_assay_values = TRUE,
    marker_col = "marker"
  )

  expect_true("R" %in% colnames(result))
  expect_true("NWT" %in% colnames(result))
  expect_true("mic" %in% colnames(result))
  expect_true(all(c("gyrA_S83F", "parC_S80I") %in% gsub("\\.\\.", ":", colnames(result))))
})

test_that("assay_by_var filters on drug via pheno_drug", {
  pheno_table <- tibble::tibble(
    id = c("s1", "s2"),
    drug = AMR::as.ab(c("Ciprofloxacin", "Ciprofloxacin")),
    mic = AMR::as.mic(c("4", "0.25")),
    pheno_clsi = AMR::as.sir(c("R", "S"))
  )

  plot_obj <- assay_by_var(
    pheno_table = pheno_table,
    pheno_drug = "Ciprofloxacin",
    measure = "mic",
    colour_by = "pheno_clsi"
  )

  expect_s3_class(plot_obj, "ggplot")
})
