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

# gtdb.mo() Tests -------------------------------------------------------

test_that("gtdb.mo works with valid GTDB species names", {
  # Test basic GTDB format with suffix
  result <- gtdb.mo("Escherichia_A coli_BC")
  expect_s3_class(result, "mo")
  expect_length(result, 1)
  
  # Test another example
  result2 <- gtdb.mo("Pseudomonas_E piscis")
  expect_s3_class(result2, "mo")
  expect_length(result2, 1)
})

test_that("gtdb.mo handles species names without GTDB suffix", {
  # Test regular species name
  result <- gtdb.mo("Escherichia coli")
  expect_s3_class(result, "mo")
  expect_length(result, 1)
})

test_that("gtdb.mo handles vector input", {
  # Test with multiple species
  species_vec <- c("Escherichia_A coli_BC", "Pseudomonas_E piscis", "Acinetobacter calcoaceticus_C")
  result <- gtdb.mo(species_vec)
  
  expect_s3_class(result, "mo")
  expect_length(result, 3)
})

test_that("gtdb.mo handles edge cases", {
  # Test with empty string
  expect_s3_class(gtdb.mo(""), "mo")
  
  # Test with NA
  result <- gtdb.mo(NA)
  expect_s3_class(result, "mo")
  expect_true(is.na(result))
})

test_that("gtdb.mo handles invalid inputs appropriately", {
  # Test with NULL should error
  expect_error(gtdb.mo(NULL))
  
  # Test with non-character input
  expect_error(gtdb.mo(123))
})

test_that("gtdb.mo properly cleans GTDB suffixes", {
  # The cleaning regex should remove _X suffixes (where X is uppercase letter(s))
  result1 <- gtdb.mo("Haemophilus_D parainfluenzae_A")
  result2 <- gtdb.mo("Haemophilus parainfluenzae")
  
  expect_s3_class(result1, "mo")
  expect_s3_class(result2, "mo")
  # Both should resolve to same or similar organism
  expect_length(result1, 1)
  expect_length(result2, 1)
})

# import_gtdb() Tests ---------------------------------------------------

test_that("import_gtdb works with data frame input", {
  # Create test data frame
  test_df <- data.frame(
    Species = c("Escherichia_A coli", "Pseudomonas_E aeruginosa"),
    SampleID = c("Sample1", "Sample2"),
    stringsAsFactors = FALSE
  )
  
  result <- import_gtdb(tbl = test_df)
  
  # Check output structure
  expect_s3_class(result, "data.frame")
  expect_true("gtdb.mo" %in% colnames(result))
  expect_true("gtdb.species" %in% colnames(result))
  expect_equal(nrow(result), 2)
  
  # Check new columns have correct classes
  expect_s3_class(result$gtdb.mo, "mo")
  expect_type(result$gtdb.species, "character")
})

test_that("import_gtdb works with custom species column", {
  # Test with different column name
  test_df <- data.frame(
    CustomSpeciesCol = c("Escherichia_A coli", "Staphylococcus aureus"),
    SampleID = c("Sample1", "Sample2"),
    stringsAsFactors = FALSE
  )
  
  result <- import_gtdb(tbl = test_df, species_column = "CustomSpeciesCol")
  
  expect_s3_class(result, "data.frame")
  expect_true("gtdb.mo" %in% colnames(result))
  expect_true("gtdb.species" %in% colnames(result))
  expect_equal(nrow(result), 2)
})

test_that("import_gtdb handles edge cases", {
  # Test with single row
  test_df <- data.frame(
    Species = "Escherichia_A coli",
    stringsAsFactors = FALSE
  )
  
  result <- import_gtdb(tbl = test_df)
  expect_equal(nrow(result), 1)
  expect_true("gtdb.mo" %in% colnames(result))
  
  # Test with empty data frame
  empty_df <- data.frame(Species = character(0))
  result_empty <- import_gtdb(tbl = empty_df)
  expect_equal(nrow(result_empty), 0)
  expect_true("gtdb.mo" %in% colnames(result_empty))
})

test_that("import_gtdb errors appropriately with invalid inputs", {
  # Test with NULL inputs
  expect_error(import_gtdb(file = NULL, tbl = NULL))
  
  # Test with missing species column
  test_df <- data.frame(NotSpecies = c("test1", "test2"))
  expect_error(import_gtdb(tbl = test_df, species_column = "Species"))
  
  # Test with non-data frame input
  expect_error(import_gtdb(tbl = "not a data frame"))
})

test_that("import_gtdb preserves original columns", {
  # Test that original columns are kept
  test_df <- data.frame(
    Species = c("Escherichia_A coli", "Pseudomonas_E aeruginosa"),
    SampleID = c("Sample1", "Sample2"),
    Coverage = c(95.5, 98.2),
    stringsAsFactors = FALSE
  )
  
  result <- import_gtdb(tbl = test_df)
  
  # Check all original columns are present
  expect_true("Species" %in% colnames(result))
  expect_true("SampleID" %in% colnames(result))
  expect_true("Coverage" %in% colnames(result))
  
  # Check new columns added
  expect_true("gtdb.mo" %in% colnames(result))
  expect_true("gtdb.species" %in% colnames(result))
  
  # Check data integrity
  expect_equal(result$SampleID, c("Sample1", "Sample2"))
  expect_equal(result$Coverage, c(95.5, 98.2))
})

test_that("import_gtdb handles various GTDB naming patterns", {
  # Test with various GTDB suffixes
  test_df <- data.frame(
    Species = c(
      "Pseudomonas_E piscis",
      "Haemophilus_D parainfluenzae_A",
      "Acinetobacter calcoaceticus_C",
      "Escherichia coli"  # No suffix
    ),
    stringsAsFactors = FALSE
  )
  
  result <- import_gtdb(tbl = test_df)
  
  expect_equal(nrow(result), 4)
  expect_s3_class(result$gtdb.mo, "mo")
  expect_length(result$gtdb.mo, 4)
  
  # All should resolve to something (even if UNKNOWN)
  expect_length(result$gtdb.species, 4)
})

test_that("import_gtdb output integrates with AMR package functions", {
  # Test that output can be used with AMR package
  test_df <- data.frame(
    Species = c("Escherichia_A coli", "Staphylococcus aureus"),
    stringsAsFactors = FALSE
  )
  
  result <- import_gtdb(tbl = test_df)
  
  # Should be able to use AMR functions on the mo column
  expect_no_error({
    AMR::mo_name(result$gtdb.mo)
    AMR::mo_genus(result$gtdb.mo)
  })
})
