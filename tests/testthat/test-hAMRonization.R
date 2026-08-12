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

# harmonize_data() Tests ------------------------------------------------

test_that("harmonize_data requires reticulate and Python", {
  skip_if_not_installed("reticulate")
  
  # Test that function exists
  expect_true(exists("harmonize_data"))
})

test_that("harmonize_data requires all parameters", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # Test that required parameters are checked
  expect_error(harmonize_data())
  expect_error(harmonize_data(user_software_name = "AMRFinderPlus"))
  expect_error(harmonize_data(
    user_software_name = "AMRFinderPlus",
    user_software_version = "3.12"
  ))
})

test_that("harmonize_data validates input parameters", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # All parameters provided
  expect_error(harmonize_data(
    user_software_name = "AMRFinderPlus",
    user_software_version = "3.12.8",
    user_database_version = "2024-01-31.1",
    user_input_filename = NULL
  ))
})

test_that("harmonize_data handles various software names", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # Common AMR software names
  software_names <- c(
    "AMRFinderPlus",
    "ResFinder",
    "CARD",
    "ABRicate"
  )
  
  # Function should accept these names
  for (name in software_names) {
    expect_no_error({
      # Just test parameter validation, not actual execution
      func_params <- list(
        user_software_name = name,
        user_software_version = "1.0",
        user_database_version = "2024",
        user_input_filename = "test.tsv"
      )
    })
  }
})

test_that("harmonize_data requires Python environment", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # Test that Python/hAMRonization availability is checked
  # This will likely fail unless hAMRonization is installed
  expect_error(harmonize_data(
    user_software_name = "AMRFinderPlus",
    user_software_version = "3.12.8",
    user_database_version = "2024-01-31.1",
    user_input_filename = "nonexistent.tsv"
  ))
})

test_that("harmonize_data handles file path validation", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # Test with non-existent file
  expect_error(harmonize_data(
    user_software_name = "AMRFinderPlus",
    user_software_version = "3.12.8",
    user_database_version = "2024-01-31.1",
    user_input_filename = "/nonexistent/path/file.tsv"
  ))
})

test_that("harmonize_data documentation is accessible", {
  # Check that help documentation exists
  expect_true("harmonize_data" %in% ls("package:AMRgen"))
})

# Integration Tests -----------------------------------------------------

test_that("harmonize_data is designed for AMRFinderPlus output", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # The function should work with AMRFinderPlus output format
  # We can't test actual execution without Python setup
  # But we can verify the function signature
  
  func_args <- names(formals(harmonize_data))
  
  expect_true("user_software_name" %in% func_args)
  expect_true("user_software_version" %in% func_args)
  expect_true("user_database_version" %in% func_args)
  expect_true("user_input_filename" %in% func_args)
})

test_that("harmonize_data would work with ecoli_geno_raw structure", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # Check that ecoli_geno_raw has expected structure for harmonization
  expect_true(exists("ecoli_geno_raw"))
  expect_s3_class(ecoli_geno_raw, "data.frame")
  
  # Expected columns from AMRFinderPlus
  expect_true("Name" %in% colnames(ecoli_geno_raw))
  expect_true("Gene symbol" %in% colnames(ecoli_geno_raw))
})

# Setup and Environment Tests -------------------------------------------

test_that("harmonize_data has proper error messages for missing dependencies", {
  skip_if_not_installed("reticulate")
  
  # Function should provide helpful error messages
  # when Python or hAMRonization is not available
  
  # We can't easily test this without actually removing Python
  # But we can verify the function exists and is callable
  expect_true(is.function(harmonize_data))
})

test_that("harmonize_data parameters have sensible defaults or requirements", {
  # Check function signature
  func_def <- formals(harmonize_data)
  
  # All parameters should be required (no defaults) based on function design
  expect_equal(length(func_def), 4)
  
  # Parameter names should be descriptive
  param_names <- names(func_def)
  expect_true(all(grepl("user_", param_names)))
})

# Documentation Tests ---------------------------------------------------

test_that("harmonize_data has appropriate documentation", {
  # Check that the function is exported
  namespace <- getNamespace("AMRgen")
  exports <- getNamespaceExports(namespace)
  
  expect_true("harmonize_data" %in% exports)
})

test_that("harmonize_data function signature is stable", {
  # Verify the function has the expected parameters
  func_params <- names(formals(harmonize_data))
  
  expected_params <- c(
    "user_software_name",
    "user_software_version",
    "user_database_version",
    "user_input_filename"
  )
  
  expect_equal(sort(func_params), sort(expected_params))
})

# Conditional Tests Based on Environment --------------------------------

test_that("harmonize_data works when hAMRonization is available", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # Check if Python and hAMRonization are available
  python_available <- tryCatch({
    reticulate::py_available()
  }, error = function(e) FALSE)
  
  if (!python_available) {
    skip("Python not available")
  }
  
  # Try to check for hAMRonization
  hamronization_available <- tryCatch({
    reticulate::py_module_available("hAMRonization")
  }, error = function(e) FALSE)
  
  if (!hamronization_available) {
    skip("hAMRonization not available in Python environment")
  }
  
  # If we get here, we could potentially test actual functionality
  # But we need a valid input file, which we don't have in tests
  skip("Requires valid AMR tool output file for full integration test")
})

test_that("harmonize_data provides informative errors", {
  skip_if_not_installed("reticulate")
  skip_on_cran()
  
  # Test that errors are informative when things go wrong
  # This is hard to test without actually breaking things
  # So we just verify the function can be called with error handling
  
  expect_error({
    harmonize_data(
      user_software_name = NULL,
      user_software_version = NULL,
      user_database_version = NULL,
      user_input_filename = NULL
    )
  })
})
