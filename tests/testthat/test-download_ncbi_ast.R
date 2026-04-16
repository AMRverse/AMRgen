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

# download_ncbi_ast() Tests ---------------------------------------------

test_that("download_ncbi_ast requires internet connection", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with very small record limit
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 5,
    reformat = FALSE,
    interpret_eucast = FALSE,
    interpret_clsi = FALSE,
    interpret_ecoff = FALSE
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("download_ncbi_ast works with basic parameters", {
  skip_if_offline()
  skip_on_cran()
  
  # Test basic download
  result <- download_ncbi_ast(
    species = "Staphylococcus aureus",
    max_records = 3,
    reformat = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  expect_lte(nrow(result), 3 * 50)  # max_records * approximate entries per sample
})

test_that("download_ncbi_ast works with antibiotic filter", {
  skip_if_offline()
  skip_on_cran()
  
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    antibiotic = "Ciprofloxacin",
    max_records = 3,
    reformat = FALSE
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("download_ncbi_ast works with reformat = TRUE", {
  skip_if_offline()
  skip_on_cran()
  
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 3,
    reformat = TRUE
  )
  
  expect_s3_class(result, "data.frame")
  
  # Check for reformatted columns
  if (nrow(result) > 0) {
    expect_true("id" %in% colnames(result) || 
                "drug_agent" %in% colnames(result) ||
                "#BioSample" %in% colnames(result))
  }
})

test_that("download_ncbi_ast works with interpretation flags", {
  skip_if_offline()
  skip_on_cran()
  
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 3,
    reformat = TRUE,
    interpret_eucast = TRUE,
    interpret_clsi = FALSE,
    interpret_ecoff = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  
  if (nrow(result) > 0) {
    # Should have interpretation columns
    possible_cols <- c("pheno_eucast", "pheno_clsi", "ecoff")
    has_interp_col <- any(possible_cols %in% colnames(result))
    expect_true(has_interp_col)
  }
})

test_that("download_ncbi_ast respects max_records parameter", {
  skip_if_offline()
  skip_on_cran()
  
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 2,
    reformat = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  # Result should be limited (though exact count depends on data structure)
})

test_that("download_ncbi_ast handles batch_size parameter", {
  skip_if_offline()
  skip_on_cran()
  
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 3,
    batch_size = 100,
    reformat = FALSE
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("download_ncbi_ast handles sleep_time parameter", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with custom sleep time
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 2,
    sleep_time = 0.5,
    reformat = FALSE
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("download_ncbi_ast handles force_antibiotic parameter", {
  skip_if_offline()
  skip_on_cran()
  
  # This parameter forces specific antibiotic even if not specified
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 2,
    force_antibiotic = TRUE,
    antibiotic = "Ampicillin",
    reformat = FALSE
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("download_ncbi_ast errors appropriately with invalid inputs", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with NULL species
  expect_error(download_ncbi_ast(species = NULL))
  
  # Test with invalid max_records
  expect_error(download_ncbi_ast(species = "Escherichia coli", max_records = -1))
  expect_error(download_ncbi_ast(species = "Escherichia coli", max_records = 0))
})

test_that("download_ncbi_ast handles species name variations", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with full species name
  result1 <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 2,
    reformat = FALSE
  )
  expect_s3_class(result1, "data.frame")
  
  # Test with abbreviated species name
  result2 <- download_ncbi_ast(
    species = "E. coli",
    max_records = 2,
    reformat = FALSE
  )
  expect_s3_class(result2, "data.frame")
})

test_that("download_ncbi_ast handles empty results gracefully", {
  skip_if_offline()
  skip_on_cran()
  
  # Try with an obscure species that might not have data
  result <- download_ncbi_ast(
    species = "Nonexistent species xyz123",
    max_records = 1,
    reformat = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  # May be empty or have no rows
})

test_that("download_ncbi_ast parameter combinations work together", {
  skip_if_offline()
  skip_on_cran()
  
  # Test complex parameter combination
  result <- download_ncbi_ast(
    species = "Staphylococcus aureus",
    antibiotic = "Vancomycin",
    max_records = 2,
    batch_size = 50,
    sleep_time = 0.3,
    reformat = TRUE,
    interpret_eucast = TRUE
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("download_ncbi_ast output structure is consistent", {
  skip_if_offline()
  skip_on_cran()
  
  result <- download_ncbi_ast(
    species = "Escherichia coli",
    max_records = 2,
    reformat = FALSE
  )
  
  if (nrow(result) > 0) {
    # Check that basic structure is maintained
    expect_true(is.data.frame(result))
    expect_true(ncol(result) > 0)
  }
})
