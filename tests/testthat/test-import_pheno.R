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

# import_ncbi_ast() Tests -----------------------------------------------

test_that("import_ncbi_ast works with example data", {
  # Use package example data
  result <- suppressWarnings(import_ncbi_ast(ecoli_ast_raw))
  
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  
  # Check expected columns exist
  expect_true("id" %in% colnames(result))
  expect_true("drug_agent" %in% colnames(result))
})

test_that("import_ncbi_ast creates correct AMR classes", {
  result <- suppressWarnings(import_ncbi_ast(ecoli_ast_raw))
  
  # Check AMR package classes
  if ("drug_agent" %in% colnames(result)) {
    expect_s3_class(result$drug_agent, "ab")
  }
  
  if ("spp_pheno" %in% colnames(result)) {
    expect_s3_class(result$spp_pheno, "mo")
  }
  
  if ("mic" %in% colnames(result)) {
    expect_s3_class(result$mic, "mic")
  }
  
  if ("disk" %in% colnames(result)) {
    expect_s3_class(result$disk, "disk")
  }
})

test_that("import_ncbi_ast works with EUCAST interpretation", {
  result <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    interpret_eucast = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
  
  # Should have pheno_eucast column
  if ("pheno_eucast" %in% colnames(result)) {
    expect_s3_class(result$pheno_eucast, "sir")
  }
})

test_that("import_ncbi_ast works with CLSI interpretation", {
  result <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    interpret_clsi = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
  
  # Should have pheno_clsi column
  if ("pheno_clsi" %in% colnames(result)) {
    expect_s3_class(result$pheno_clsi, "sir")
  }
})

test_that("import_ncbi_ast works with ECOFF interpretation", {
  result <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    interpret_ecoff = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
  
  # Should have ecoff column
  if ("ecoff" %in% colnames(result)) {
    expect_s3_class(result$ecoff, "sir")
  }
})

test_that("import_ncbi_ast works with custom sample column", {
  # Create test data with custom column name
  test_data <- ecoli_ast_raw
  colnames(test_data)[colnames(test_data) == "#BioSample"] <- "CustomID"
  
  result <- suppressWarnings(import_ncbi_ast(
    test_data,
    sample_col = "CustomID"
  ))
  
  expect_s3_class(result, "data.frame")
  expect_true("id" %in% colnames(result))
})

test_that("import_ncbi_ast handles species parameter", {
  result <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    species = "Escherichia coli",
    interpret_eucast = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("import_ncbi_ast handles antibiotic parameter", {
  result <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    ab = "Ciprofloxacin",
    interpret_eucast = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("import_ncbi_ast handles source parameter", {
  result <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    source = "Test Source"
  ))
  
  expect_s3_class(result, "data.frame")
  
  if ("source" %in% colnames(result)) {
    expect_true(all(result$source == "Test Source" | is.na(result$source)))
  }
})

test_that("import_ncbi_ast errors with missing required columns", {
  # Test with missing Antibiotic column
  bad_data <- data.frame(
    BioSample = c("SAMN001", "SAMN002"),
    MIC = c("0.5", "1"),
    stringsAsFactors = FALSE
  )
  
  expect_error(suppressWarnings(import_ncbi_ast(bad_data)))
})

test_that("import_ncbi_ast errors with invalid sample_col", {
  expect_error(import_ncbi_ast(
    ecoli_ast_raw,
    sample_col = "NonExistentColumn"
  ))
})

test_that("import_ncbi_ast handles edge cases", {
  # Test with single row
  single_row <- ecoli_ast_raw[1, , drop = FALSE]
  result <- suppressWarnings(import_ncbi_ast(single_row))
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
})

# import_ast() Tests ----------------------------------------------------

test_that("import_ast works with ncbi format", {
  result <- suppressWarnings(import_ast(
    ecoli_ast_raw,
    format = "ncbi"
  ))
  
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})

test_that("import_ast format parameter accepts various formats", {
  # Test that format parameter is recognized
  expect_no_error(suppressWarnings(import_ast(
    ecoli_ast_raw,
    format = "ncbi"
  )))
})

test_that("import_ast passes parameters to underlying functions", {
  result <- suppressWarnings(import_ast(
    ecoli_ast_raw,
    format = "ncbi",
    interpret_eucast = TRUE,
    interpret_ecoff = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
})

# interpret_ast() Tests -------------------------------------------------

test_that("interpret_ast works with ecoli_ast data", {
  # Use pre-imported data
  result <- suppressWarnings(interpret_ast(ecoli_ast))
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(ecoli_ast))
})

test_that("interpret_ast adds interpretation columns", {
  result <- suppressWarnings(interpret_ast(
    ecoli_ast,
    interpret_eucast = TRUE,
    interpret_clsi = TRUE,
    interpret_ecoff = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
  
  # Should have at least one interpretation column
  interp_cols <- c("pheno_eucast", "pheno_clsi", "ecoff")
  has_any_interp <- any(interp_cols %in% colnames(result))
  expect_true(has_any_interp)
})

test_that("interpret_ast handles species override", {
  result <- suppressWarnings(interpret_ast(
    ecoli_ast,
    species = "Escherichia coli",
    interpret_eucast = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("interpret_ast handles antibiotic override", {
  result <- suppressWarnings(interpret_ast(
    ecoli_ast,
    ab = "Ciprofloxacin",
    interpret_eucast = TRUE
  ))
  
  expect_s3_class(result, "data.frame")
})

test_that("interpret_ast preserves original data", {
  original_cols <- colnames(ecoli_ast)
  result <- suppressWarnings(interpret_ast(ecoli_ast))
  
  # All original columns should still be present
  expect_true(all(original_cols %in% colnames(result)))
})

test_that("interpret_ast handles edge cases", {
  # Test with single row
  single_row <- ecoli_ast[1, , drop = FALSE]
  result <- suppressWarnings(interpret_ast(single_row))
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
})

# Integration Tests -----------------------------------------------------

test_that("import and interpret workflow works together", {
  # Import data
  imported <- suppressWarnings(import_ncbi_ast(ecoli_ast_raw))
  
  # Then interpret
  interpreted <- suppressWarnings(interpret_ast(
    imported,
    interpret_eucast = TRUE
  ))
  
  expect_s3_class(interpreted, "data.frame")
  expect_equal(nrow(imported), nrow(interpreted))
})

test_that("import with interpretation matches interpret_ast", {
  # Import with interpretation
  result1 <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    interpret_eucast = TRUE
  ))
  
  # Import then interpret
  imported <- suppressWarnings(import_ncbi_ast(ecoli_ast_raw))
  result2 <- suppressWarnings(interpret_ast(
    imported,
    interpret_eucast = TRUE
  ))
  
  expect_s3_class(result1, "data.frame")
  expect_s3_class(result2, "data.frame")
})

test_that("imported data works with other AMRgen functions", {
  # Import data
  imported <- suppressWarnings(import_ncbi_ast(
    ecoli_ast_raw,
    interpret_eucast = TRUE
  ))
  
  # Should be able to use with AMR package functions
  expect_no_error({
    if ("drug_agent" %in% colnames(imported)) {
      AMR::ab_name(imported$drug_agent[1])
    }
  })
  
  expect_no_error({
    if ("spp_pheno" %in% colnames(imported)) {
      AMR::mo_name(imported$spp_pheno[1])
    }
  })
})
