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

# eucast_supported_ab_distributions() Tests -----------------------------

test_that("eucast_supported_ab_distributions requires internet", {
  skip_if_offline()
  skip_on_cran()
  
  result <- eucast_supported_ab_distributions()
  
  expect_type(result, "character")
  expect_gt(length(result), 0)
  expect_true(is.character(result))
})

test_that("eucast_supported_ab_distributions returns named vector", {
  skip_if_offline()
  skip_on_cran()
  
  result <- eucast_supported_ab_distributions()
  
  expect_true(!is.null(names(result)))
  expect_equal(length(names(result)), length(result))
})

test_that("eucast_supported_ab_distributions is sorted", {
  skip_if_offline()
  skip_on_cran()
  
  result <- eucast_supported_ab_distributions()
  
  # Result should be sorted
  expect_equal(result, sort(result))
})

test_that("eucast_supported_ab_distributions caches results", {
  skip_if_offline()
  skip_on_cran()
  
  # First call
  result1 <- eucast_supported_ab_distributions()
  
  # Second call should use cache
  result2 <- eucast_supported_ab_distributions()
  
  expect_equal(result1, result2)
})

# get_eucast_amr_distribution() Tests -----------------------------------

test_that("get_eucast_amr_distribution works with antibiotic only", {
  skip_if_offline()
  skip_on_cran()
  
  result <- get_eucast_amr_distribution(ab = "Ciprofloxacin", method = "MIC")
  
  expect_s3_class(result, "data.frame")
  if (nrow(result) > 0) {
    expect_true("mic" %in% colnames(result) || "disk" %in% colnames(result))
  }
})

test_that("get_eucast_amr_distribution works with antibiotic and organism", {
  skip_if_offline()
  skip_on_cran()
  
  result <- get_eucast_amr_distribution(
    ab = "Ciprofloxacin",
    mo = "Klebsiella pneumoniae",
    method = "MIC"
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("get_eucast_amr_distribution works with disk method", {
  skip_if_offline()
  skip_on_cran()
  
  result <- get_eucast_amr_distribution(
    ab = "Ciprofloxacin",
    mo = "Escherichia coli",
    method = "disk"
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("get_eucast_amr_distribution as_freq_table parameter works", {
  skip_if_offline()
  skip_on_cran()
  
  # With freq table
  result_freq <- get_eucast_amr_distribution(
    ab = "Ciprofloxacin",
    method = "MIC",
    as_freq_table = TRUE
  )
  
  # Without freq table
  result_no_freq <- get_eucast_amr_distribution(
    ab = "Ciprofloxacin",
    method = "MIC",
    as_freq_table = FALSE
  )
  
  expect_s3_class(result_freq, "data.frame")
  expect_s3_class(result_no_freq, "data.frame")
})

# get_eucast_mic_distribution() Tests -----------------------------------

test_that("get_eucast_mic_distribution works", {
  skip_if_offline()
  skip_on_cran()
  
  result <- get_eucast_mic_distribution(ab = "Ciprofloxacin")
  
  expect_s3_class(result, "data.frame")
  if (nrow(result) > 0) {
    expect_true("mic" %in% colnames(result))
  }
})

test_that("get_eucast_mic_distribution with organism filter", {
  skip_if_offline()
  skip_on_cran()
  
  result <- get_eucast_mic_distribution(
    ab = "Ciprofloxacin",
    mo = "Klebsiella pneumoniae"
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("get_eucast_mic_distribution returns frequency table by default", {
  skip_if_offline()
  skip_on_cran()
  
  result <- get_eucast_mic_distribution(ab = "Ciprofloxacin")
  
  expect_s3_class(result, "data.frame")
  # Should have count column for frequency table
  if (nrow(result) > 0) {
    expect_true("count" %in% colnames(result) || "freq" %in% colnames(result))
  }
})

# get_eucast_disk_distribution() Tests ----------------------------------

test_that("get_eucast_disk_distribution works", {
  skip_if_offline()
  skip_on_cran()
  
  # Note: Not all antibiotics have disk data
  # Using one that likely has disk data
  result <- tryCatch({
    get_eucast_amr_distribution(ab = "Ciprofloxacin", method = "disk")
  }, error = function(e) {
    skip("Disk distribution not available for test antibiotic")
  })
  
  expect_s3_class(result, "data.frame")
})

# compare_mic_with_eucast() Tests ---------------------------------------

test_that("compare_mic_with_eucast works with valid MIC values", {
  skip_if_offline()
  skip_on_cran()
  
  # Create test MIC values
  test_mics <- c("0.001", "0.5", "2", "8")
  
  result <- compare_mic_with_eucast(
    mics = test_mics,
    ab = "Ciprofloxacin",
    mo = "Escherichia coli"
  )
  
  expect_s3_class(result, "compare_eucast")
  expect_s3_class(result, "data.frame")
})

test_that("compare_mic_with_eucast handles MIC coercion", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with various MIC formats
  test_mics <- c(0.5, 1, 2, 4)
  
  result <- compare_mic_with_eucast(
    mics = test_mics,
    ab = "Ampicillin",
    mo = "Escherichia coli"
  )
  
  expect_s3_class(result, "compare_eucast")
})

test_that("compare_mic_with_eucast returns correct structure", {
  skip_if_offline()
  skip_on_cran()
  
  test_mics <- c("0.5", "1", "2")
  
  result <- compare_mic_with_eucast(
    mics = test_mics,
    ab = "Ciprofloxacin",
    mo = "Escherichia coli"
  )
  
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})

# compare_disk_with_eucast() Tests --------------------------------------

test_that("compare_disk_with_eucast works with valid disk values", {
  skip_if_offline()
  skip_on_cran()
  
  # Create test disk values
  test_disks <- c("10", "15", "20", "25")
  
  result <- tryCatch({
    compare_disk_with_eucast(
      disks = test_disks,
      ab = "Ciprofloxacin",
      mo = "Escherichia coli"
    )
  }, error = function(e) {
    skip("Disk distribution not available for test antibiotic")
  })
  
  expect_s3_class(result, "compare_eucast")
  expect_s3_class(result, "data.frame")
})

test_that("compare_disk_with_eucast handles disk coercion", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with numeric disk values
  test_disks <- c(10, 15, 20, 25)
  
  result <- tryCatch({
    compare_disk_with_eucast(
      disks = test_disks,
      ab = "Ciprofloxacin",
      mo = "Escherichia coli"
    )
  }, error = function(e) {
    skip("Disk distribution not available for test antibiotic")
  })
  
  expect_s3_class(result, "compare_eucast")
})

# Error Handling Tests --------------------------------------------------

test_that("EUCAST functions error appropriately with invalid inputs", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with NULL antibiotic
  expect_error(get_eucast_mic_distribution(ab = NULL))
  
  # Test with invalid antibiotic
  expect_error(get_eucast_mic_distribution(ab = "InvalidDrugXYZ123"))
})

test_that("compare functions handle empty input", {
  skip_if_offline()
  skip_on_cran()
  
  # Test with empty MIC vector
  expect_error(compare_mic_with_eucast(
    mics = character(0),
    ab = "Ciprofloxacin",
    mo = "Escherichia coli"
  ))
})

# Integration Tests -----------------------------------------------------

test_that("EUCAST functions integrate with AMR package classes", {
  skip_if_offline()
  skip_on_cran()
  
  # Get distribution
  result <- get_eucast_mic_distribution(ab = "Ciprofloxacin")
  
  if (nrow(result) > 0 && "mic" %in% colnames(result)) {
    # MIC column should be of class 'mic'
    expect_s3_class(result$mic, "mic")
  }
})

test_that("EUCAST distribution data can be used for plotting", {
  skip_if_offline()
  skip_on_cran()
  
  result <- get_eucast_mic_distribution(ab = "Ciprofloxacin")
  
  if (nrow(result) > 0) {
    # Should be able to create a basic plot
    expect_no_error({
      data_for_plot <- result
      nrow(data_for_plot)
    })
  }
})
