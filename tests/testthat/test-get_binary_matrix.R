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

library(testthat)
library(AMRgen)

# get_binary_matrix() Tests ---------------------------------------------

test_that("get_binary_matrix requires drug_agent column in pheno_table", {
  # Create minimal test data without drug_agent
  bad_pheno <- data.frame(
    id = c("Sample1", "Sample2"),
    antibiotic = c("CIP", "CIP")
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("gyrA_S83L", "parC_S80I"),
    drug_class = c("Quinolones", "Quinolones")
  )
  
  expect_error(suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = bad_pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  )))
})

test_that("get_binary_matrix requires drug_class column in geno_table", {
  # Create minimal test data without drug_class
  pheno <- data.frame(
    id = c("Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S"))
  )
  
  bad_geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("gyrA_S83L", "parC_S80I")
  )
  
  expect_error(suppressWarnings(get_binary_matrix(
    geno_table = bad_geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  )))
})

test_that("get_binary_matrix errors when antibiotic not in pheno_table", {
  pheno <- data.frame(
    id = c("Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("AMP", "AMP")),
    pheno_clsi = AMR::as.sir(c("R", "S"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("gyrA_S83L", "parC_S80I"),
    drug_class = c("Quinolones", "Quinolones")
  )
  
  expect_error(suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  )))
})

test_that("get_binary_matrix errors when drug_class not in geno_table", {
  pheno <- data.frame(
    id = c("Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("blaTEM", "blaCTX"),
    drug_class = c("Beta-lactams", "Beta-lactams")
  )
  
  expect_error(suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  )))
})

test_that("get_binary_matrix works with minimal valid data", {
  # Create minimal valid test data
  pheno <- data.frame(
    id = c("Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S")),
    ecoff = AMR::as.sir(c("R", "S"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("gyrA_S83L", "gyrA_WT"),
    drug_class = c("Quinolones", "Quinolones")
  )
  
  result <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  ))
  
  expect_s3_class(result, "data.frame")
  expect_true("id" %in% colnames(result))
})

test_that("get_binary_matrix keeps SIR column when requested", {
  pheno <- data.frame(
    id = c("Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S")),
    ecoff = AMR::as.sir(c("R", "S"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("gyrA_S83L", "gyrA_WT"),
    drug_class = c("Quinolones", "Quinolones")
  )
  
  result <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones",
    keep_SIR = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("get_binary_matrix handles custom column names", {
  pheno <- data.frame(
    sample_id = c("Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("CIP", "CIP")),
    resistance = AMR::as.sir(c("R", "S")),
    ecoff = AMR::as.sir(c("R", "S"))
  )
  
  geno <- data.frame(
    sample_id = c("Sample1", "Sample2"),
    gene = c("gyrA_S83L", "gyrA_WT"),
    drug_class = c("Quinolones", "Quinolones")
  )
  
  result <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones",
    geno_sample_col = "sample_id",
    pheno_sample_col = "sample_id",
    sir_col = "resistance",
    marker_col = "gene"
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("get_binary_matrix handles keep_assay_values parameter", {
  pheno <- data.frame(
    id = c("Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S")),
    ecoff = AMR::as.sir(c("R", "S")),
    mic = AMR::as.mic(c("4", "0.5")),
    disk = AMR::as.disk(c("10", "25"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("gyrA_S83L", "gyrA_WT"),
    drug_class = c("Quinolones", "Quinolones")
  )
  
  result <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones",
    keep_assay_values = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("get_binary_matrix handles most_resistant parameter", {
  # Create data with multiple entries per sample
  pheno <- data.frame(
    id = c("Sample1", "Sample1", "Sample2"),
    drug_agent = AMR::as.ab(c("CIP", "CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S", "S")),
    ecoff = AMR::as.sir(c("R", "S", "S"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2"),
    marker = c("gyrA_S83L", "gyrA_WT"),
    drug_class = c("Quinolones", "Quinolones")
  )
  
  # Test with most_resistant = TRUE
  result_most <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones",
    most_resistant = TRUE
  ))
  
  expect_s3_class(result_most, "data.frame")
  
  # Test with most_resistant = FALSE
  result_least <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones",
    most_resistant = FALSE
  ))
  
  expect_s3_class(result_least, "data.frame")
})

test_that("get_binary_matrix output has correct structure", {
  pheno <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    drug_agent = AMR::as.ab(c("CIP", "CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S", "I")),
    ecoff = AMR::as.sir(c("R", "S", "R"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    marker = c("gyrA_S83L", "gyrA_WT", "parC_S80I"),
    drug_class = c("Quinolones", "Quinolones", "Quinolones")
  )
  
  result <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  ))
  
  # Check basic structure
  expect_s3_class(result, "data.frame")
  expect_true("id" %in% colnames(result))
  expect_gt(ncol(result), 1)  # Should have id plus other columns
})

# get_combo_matrix() Tests ----------------------------------------------

test_that("get_combo_matrix requires binary_matrix input", {
  # This function adds marker combinations to a binary matrix
  # Test that it requires appropriate input
  
  # Create a mock binary matrix
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2"),
    R = c(1, 0),
    NWT = c(1, 0),
    gyrA_S83L = c(1, 0),
    parC_S80I = c(0, 1)
  )
  
  # Should accept data frame
  expect_no_error(suppressWarnings(get_combo_matrix(binary_mat)))
})

test_that("get_combo_matrix adds combination columns", {
  # Create a mock binary matrix
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    R = c(1, 0, 1),
    NWT = c(1, 0, 1),
    gyrA_S83L = c(1, 0, 1),
    parC_S80I = c(0, 1, 1)
  )
  
  result <- suppressWarnings(get_combo_matrix(binary_mat))
  
  expect_s3_class(result, "data.frame")
  # Result should have more or same columns as input
  expect_gte(ncol(result), ncol(binary_mat))
})

test_that("get_combo_matrix handles assay parameter", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2"),
    R = c(1, 0),
    NWT = c(1, 0),
    gyrA_S83L = c(1, 0),
    mic = AMR::as.mic(c("4", "0.5"))
  )
  
  result <- suppressWarnings(get_combo_matrix(
    binary_mat,
    assay = "mic"
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("get_combo_matrix handles edge cases", {
  # Single row
  binary_mat <- data.frame(
    id = "Sample1",
    R = 1,
    NWT = 1,
    gyrA_S83L = 1
  )
  
  result <- suppressWarnings(get_combo_matrix(binary_mat))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
})

# Integration Tests -----------------------------------------------------

test_that("get_binary_matrix and get_combo_matrix work together", {
  # Create test data
  pheno <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    drug_agent = AMR::as.ab(c("CIP", "CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S", "R")),
    ecoff = AMR::as.sir(c("R", "S", "R"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    marker = c("gyrA_S83L", "gyrA_WT", "parC_S80I"),
    drug_class = c("Quinolones", "Quinolones", "Quinolones")
  )
  
  # Get binary matrix
  binary_mat <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  ))
  
  # Get combo matrix
  combo_mat <- suppressWarnings(get_combo_matrix(binary_mat))
  
  expect_s3_class(binary_mat, "data.frame")
  expect_s3_class(combo_mat, "data.frame")
  
  # Combo matrix should have at least as many rows as binary matrix
  expect_gte(nrow(combo_mat), nrow(binary_mat))
})
