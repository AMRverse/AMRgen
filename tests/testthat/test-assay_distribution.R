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

# assay_by_var() Tests --------------------------------------------------

test_that("assay_by_var returns a ggplot object", {
  # Use package example data
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var works with disk measure", {
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "disk"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var handles colour_by parameter", {
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic",
    colour_by = "pheno_clsi"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var handles facet_var parameter", {
  # Add a faceting variable to test data
  test_data <- ecoli_ast
  test_data$source <- sample(c("Source1", "Source2"), nrow(test_data), replace = TRUE)
  
  result <- assay_by_var(
    pheno_table = test_data,
    antibiotic = "Ciprofloxacin",
    measure = "mic",
    facet_var = "source"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var handles breakpoint parameters", {
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic",
    species = "Escherichia coli",
    guideline = "EUCAST 2024",
    bp_site = NULL
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var handles custom breakpoints", {
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic",
    bp_S = 0.5,
    bp_R = 1
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var handles ECOFF parameter", {
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic",
    bp_ecoff = 0.25
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var errors with missing required parameters", {
  # Test without pheno_table
  expect_error(assay_by_var(antibiotic = "Ciprofloxacin"))
  
  # Test without antibiotic
  expect_error(assay_by_var(pheno_table = ecoli_ast))
})

test_that("assay_by_var handles edge cases", {
  # Filter to small subset
  small_data <- ecoli_ast[1:5, ]
  
  result <- assay_by_var(
    pheno_table = small_data,
    antibiotic = "Ciprofloxacin",
    measure = "mic"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var handles custom bar colors", {
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic",
    bar_cols = c("red", "blue", "green")
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var works with different antibiotics", {
  # Test with different antibiotic
  result <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ampicillin",
    measure = "mic"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var handles missing data appropriately", {
  # Create data with some missing values
  test_data <- ecoli_ast
  test_data$mic[1:3] <- NA
  
  result <- assay_by_var(
    pheno_table = test_data,
    antibiotic = "Ciprofloxacin",
    measure = "mic"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("assay_by_var plot can be modified", {
  # Test that returned plot can be further customized
  base_plot <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic"
  )
  
  # Should be able to add ggplot2 layers
  expect_no_error({
    modified_plot <- base_plot + ggplot2::labs(title = "Modified Plot")
  })
})

test_that("assay_by_var handles various measure types", {
  # MIC
  result_mic <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic"
  )
  expect_s3_class(result_mic, "ggplot")
  
  # Disk
  result_disk <- assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "disk"
  )
  expect_s3_class(result_disk, "ggplot")
})

test_that("assay_by_var integrates with breakpoint functions", {
  # Should work with getBreakpoints output
  result <- suppressMessages(assay_by_var(
    pheno_table = ecoli_ast,
    antibiotic = "Ciprofloxacin",
    measure = "mic",
    species = "Escherichia coli",
    guideline = "EUCAST 2024"
  ))
  
  expect_s3_class(result, "ggplot")
})
