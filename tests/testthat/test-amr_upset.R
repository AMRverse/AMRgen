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

# combo_stats() Tests ---------------------------------------------------

test_that("combo_stats requires binary_matrix input", {
  # Test that it requires appropriate input
  expect_error(combo_stats(NULL))
  expect_error(combo_stats("not a data frame"))
})

test_that("combo_stats works with minimal binary matrix", {
  # Create minimal test binary matrix
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    R = c(1, 0, 1),
    NWT = c(1, 0, 1),
    gyrA_S83L = c(1, 0, 1),
    parC_S80I = c(0, 1, 1)
  )
  
  result <- suppressWarnings(combo_stats(binary_mat))
  
  # Should return a list
  expect_type(result, "list")
})

test_that("combo_stats returns expected components", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    R = c(1, 0, 1, 1),
    NWT = c(1, 0, 1, 1),
    gyrA_S83L = c(1, 0, 1, 0),
    parC_S80I = c(0, 1, 1, 1)
  )
  
  result <- suppressWarnings(combo_stats(binary_mat, min_set_size = 1))
  
  expect_type(result, "list")
  # Should have various plot components
  expect_true(length(result) > 0)
})

test_that("combo_stats handles min_set_size parameter", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    R = c(1, 0, 1, 1),
    NWT = c(1, 0, 1, 1),
    gyrA_S83L = c(1, 0, 1, 0),
    parC_S80I = c(0, 1, 1, 1)
  )
  
  # Test with different min_set_size values
  result1 <- suppressWarnings(combo_stats(binary_mat, min_set_size = 1))
  result2 <- suppressWarnings(combo_stats(binary_mat, min_set_size = 2))
  
  expect_type(result1, "list")
  expect_type(result2, "list")
})

test_that("combo_stats handles order parameter", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    R = c(1, 0, 1),
    NWT = c(1, 0, 1),
    gyrA_S83L = c(1, 0, 1),
    parC_S80I = c(0, 1, 1)
  )
  
  result <- suppressWarnings(combo_stats(
    binary_mat,
    order = "freq"
  ))
  
  expect_type(result, "list")
})

test_that("combo_stats handles assay parameter", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3"),
    R = c(1, 0, 1),
    NWT = c(1, 0, 1),
    gyrA_S83L = c(1, 0, 1),
    parC_S80I = c(0, 1, 1),
    mic = AMR::as.mic(c("4", "0.5", "2"))
  )
  
  result <- suppressWarnings(combo_stats(
    binary_mat,
    assay = "mic"
  ))
  
  expect_type(result, "list")
})

test_that("combo_stats handles edge cases", {
  # Single sample
  single_sample <- data.frame(
    id = "Sample1",
    R = 1,
    NWT = 1,
    gyrA_S83L = 1,
    parC_S80I = 0
  )
  
  result <- suppressWarnings(combo_stats(single_sample, min_set_size = 1))
  expect_type(result, "list")
})

# amr_upset() Tests -----------------------------------------------------

test_that("amr_upset works with binary_matrix input", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    R = c(1, 0, 1, 1),
    NWT = c(1, 0, 1, 1),
    gyrA_S83L = c(1, 0, 1, 0),
    parC_S80I = c(0, 1, 1, 1)
  )
  
  result <- suppressWarnings(amr_upset(binary_mat, min_set_size = 1))
  
  # Should return a plot object
  expect_true(inherits(result, "ggplot") || inherits(result, "patchwork"))
})

test_that("amr_upset handles min_set_size parameter", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    R = c(1, 0, 1, 1),
    NWT = c(1, 0, 1, 1),
    gyrA_S83L = c(1, 0, 1, 0),
    parC_S80I = c(0, 1, 1, 1)
  )
  
  result <- suppressWarnings(amr_upset(
    binary_mat,
    min_set_size = 2
  ))
  
  expect_true(inherits(result, "ggplot") || inherits(result, "patchwork"))
})

test_that("amr_upset handles order parameter", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    R = c(1, 0, 1, 1),
    NWT = c(1, 0, 1, 1),
    gyrA_S83L = c(1, 0, 1, 0),
    parC_S80I = c(0, 1, 1, 1)
  )
  
  result <- suppressWarnings(amr_upset(
    binary_mat,
    order = "freq",
    min_set_size = 1
  ))
  
  expect_true(inherits(result, "ggplot") || inherits(result, "patchwork"))
})

test_that("amr_upset handles printing parameters", {
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    R = c(1, 0, 1, 1),
    NWT = c(1, 0, 1, 1),
    gyrA_S83L = c(1, 0, 1, 0),
    parC_S80I = c(0, 1, 1, 1)
  )
  
  result <- suppressWarnings(amr_upset(
    binary_mat,
    min_set_size = 1,
    print_summary = FALSE,
    print_upset = FALSE
  ))
  
  # Should still return something
  expect_true(!is.null(result))
})

# Integration Tests -----------------------------------------------------

test_that("combo_stats and amr_upset work together", {
  # Create test binary matrix
  binary_mat <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    R = c(1, 0, 1, 1),
    NWT = c(1, 0, 1, 1),
    gyrA_S83L = c(1, 0, 1, 0),
    parC_S80I = c(0, 1, 1, 1)
  )
  
  # Get combo stats
  stats <- suppressWarnings(combo_stats(binary_mat, min_set_size = 1))
  
  # Get upset plot
  upset_plot <- suppressWarnings(amr_upset(binary_mat, min_set_size = 1))
  
  expect_type(stats, "list")
  expect_true(inherits(upset_plot, "ggplot") || inherits(upset_plot, "patchwork"))
})

test_that("upset functions work with get_binary_matrix output", {
  # Create test data
  pheno <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    drug_agent = AMR::as.ab(c("CIP", "CIP", "CIP", "CIP")),
    pheno_clsi = AMR::as.sir(c("R", "S", "R", "R")),
    ecoff = AMR::as.sir(c("R", "S", "R", "R"))
  )
  
  geno <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4"),
    marker = c("gyrA_S83L", "gyrA_WT", "parC_S80I", "gyrA_S83L"),
    drug_class = c("Quinolones", "Quinolones", "Quinolones", "Quinolones")
  )
  
  # Get binary matrix
  binary_mat <- suppressWarnings(get_binary_matrix(
    geno_table = geno,
    pheno_table = pheno,
    antibiotic = "Ciprofloxacin",
    drug_class_list = "Quinolones"
  ))
  
  # Should be able to use with upset functions
  stats <- suppressWarnings(combo_stats(binary_mat, min_set_size = 1))
  upset_plot <- suppressWarnings(amr_upset(binary_mat, min_set_size = 1))
  
  expect_type(stats, "list")
  expect_true(inherits(upset_plot, "ggplot") || 
              inherits(upset_plot, "patchwork") ||
              !is.null(upset_plot))
})

test_that("upset functions handle various binary matrix structures", {
  # Test with different numbers of markers
  binary_mat_small <- data.frame(
    id = c("Sample1", "Sample2"),
    R = c(1, 0),
    NWT = c(1, 0),
    marker1 = c(1, 0)
  )
  
  binary_mat_large <- data.frame(
    id = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5"),
    R = c(1, 0, 1, 1, 0),
    NWT = c(1, 0, 1, 1, 0),
    marker1 = c(1, 0, 1, 0, 0),
    marker2 = c(0, 1, 1, 1, 0),
    marker3 = c(1, 1, 0, 1, 1)
  )
  
  result_small <- suppressWarnings(combo_stats(binary_mat_small, min_set_size = 1))
  result_large <- suppressWarnings(combo_stats(binary_mat_large, min_set_size = 1))
  
  expect_type(result_small, "list")
  expect_type(result_large, "list")
})
