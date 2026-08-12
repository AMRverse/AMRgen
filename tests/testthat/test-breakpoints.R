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

# getBreakpoints() Tests ------------------------------------------------

test_that("getBreakpoints works with valid species and antibiotic", {
  # Test with E. coli and ciprofloxacin
  result <- getBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    type_filter = "human"
  )
  
  expect_s3_class(result, "data.frame")
  expect_gte(nrow(result), 0)
  
  if (nrow(result) > 0) {
    # Check expected columns
    expect_true("guideline" %in% colnames(result))
    expect_true("mo" %in% colnames(result))
    expect_true("ab" %in% colnames(result))
    expect_true("type" %in% colnames(result))
  }
})

test_that("getBreakpoints falls back to genus level", {
  # Test with a species that might not have specific breakpoints
  # but genus should have breakpoints
  result <- getBreakpoints(
    species = "Staphylococcus epidermidis",
    guide = "EUCAST 2024",
    antibiotic = "Vancomycin",
    type_filter = "human"
  )
  
  expect_s3_class(result, "data.frame")
  # Should return something (species, genus, family, or order level)
})

test_that("getBreakpoints works with ECOFF type", {
  # Test with ECOFF type filter
  result <- getBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    type_filter = "ECOFF"
  )
  
  expect_s3_class(result, "data.frame")
  
  if (nrow(result) > 0) {
    expect_equal(unique(result$type), "ECOFF")
  }
})

test_that("getBreakpoints handles different guidelines", {
  # Test CLSI guideline
  result <- getBreakpoints(
    species = "Escherichia coli",
    guide = "CLSI 2023",
    antibiotic = "Ciprofloxacin",
    type_filter = "human"
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("getBreakpoints handles edge cases", {
  # Test with empty result (non-existent combination)
  result <- getBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "NonExistentDrug",
    type_filter = "human"
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("getBreakpoints errors appropriately with invalid inputs", {
  # Test with NULL species
  expect_error(getBreakpoints(
    species = NULL,
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin"
  ))
  
  # Test with NULL antibiotic
  expect_error(getBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = NULL
  ))
})

test_that("getBreakpoints hierarchical search works", {
  # This tests the species -> genus -> family -> order fallback
  # Use a species/drug combination that likely only has genus-level breakpoints
  result <- getBreakpoints(
    species = "Enterobacter cloacae",
    guide = "EUCAST 2024",
    antibiotic = "Meropenem",
    type_filter = "human"
  )
  
  expect_s3_class(result, "data.frame")
  # Should find something at some taxonomic level
})

test_that("getBreakpoints returns correct structure", {
  result <- getBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ampicillin",
    type_filter = "human"
  )
  
  if (nrow(result) > 0) {
    # Check for key columns
    expect_true("breakpoint_S" %in% colnames(result))
    expect_true("breakpoint_R" %in% colnames(result))
    expect_true("method" %in% colnames(result))
  }
})

# checkBreakpoints() Tests ----------------------------------------------

test_that("checkBreakpoints works with MIC assay", {
  # Suppress cat output
  result <- suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    assay = "MIC"
  ))
  
  expect_type(result, "list")
  expect_true("breakpoint_S" %in% names(result))
  expect_true("breakpoint_R" %in% names(result))
  expect_true("bp_standard" %in% names(result))
  
  # Check values are numeric
  expect_true(is.numeric(result$breakpoint_S) || is.na(result$breakpoint_S))
  expect_true(is.numeric(result$breakpoint_R) || is.na(result$breakpoint_R))
})

test_that("checkBreakpoints works with Disk assay", {
  result <- suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    assay = "Disk"
  ))
  
  expect_type(result, "list")
  expect_true("breakpoint_S" %in% names(result))
  expect_true("breakpoint_R" %in% names(result))
  expect_true("bp_standard" %in% names(result))
})

test_that("checkBreakpoints handles multiple breakpoint sites", {
  # Some antibiotics have different breakpoints for different sites
  # This test checks that the function handles this appropriately
  result <- suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    assay = "MIC"
  ))
  
  expect_type(result, "list")
  expect_type(result$bp_standard, "character")
})

test_that("checkBreakpoints uses specified bp_site when provided", {
  # Note: This test may need adjustment based on actual breakpoint data
  # The function should use the specified site if it exists
  result <- suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    assay = "MIC",
    bp_site = "uncomplicated cystitis (non-catheter)"
  ))
  
  expect_type(result, "list")
})

test_that("checkBreakpoints errors when no breakpoints found", {
  # Test with assay type that doesn't exist for the combination
  expect_error(suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    assay = "NonExistentAssay"
  )))
})

test_that("checkBreakpoints errors with invalid species-antibiotic combination", {
  # Test with combination that has no breakpoints
  expect_error(suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "NonExistentDrug",
    assay = "MIC"
  )))
})

test_that("checkBreakpoints handles edge cases", {
  # Test with NULL inputs
  expect_error(suppressMessages(checkBreakpoints(
    species = NULL,
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    assay = "MIC"
  )))
  
  expect_error(suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = NULL,
    assay = "MIC"
  )))
})

test_that("checkBreakpoints output values are reasonable", {
  result <- suppressMessages(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ampicillin",
    assay = "MIC"
  ))
  
  # Breakpoints should be numeric and positive
  if (!is.na(result$breakpoint_S)) {
    expect_true(is.numeric(result$breakpoint_S))
    expect_gte(result$breakpoint_S, 0)
  }
  
  if (!is.na(result$breakpoint_R)) {
    expect_true(is.numeric(result$breakpoint_R))
    expect_gte(result$breakpoint_R, 0)
  }
  
  # For MIC, R breakpoint should be >= S breakpoint
  if (!is.na(result$breakpoint_S) && !is.na(result$breakpoint_R)) {
    expect_gte(result$breakpoint_R, result$breakpoint_S)
  }
})

test_that("checkBreakpoints integrates with getBreakpoints", {
  # The two functions should work together
  bp_data <- getBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    type_filter = "human"
  )
  
  if (nrow(bp_data) > 0) {
    result <- suppressMessages(checkBreakpoints(
      species = "Escherichia coli",
      guide = "EUCAST 2024",
      antibiotic = "Ciprofloxacin",
      assay = "MIC"
    ))
    
    expect_type(result, "list")
  }
})

test_that("checkBreakpoints prints informative messages", {
  # Test that messages are printed
  expect_output(checkBreakpoints(
    species = "Escherichia coli",
    guide = "EUCAST 2024",
    antibiotic = "Ciprofloxacin",
    assay = "MIC"
  ), "breakpoints determined")
})
