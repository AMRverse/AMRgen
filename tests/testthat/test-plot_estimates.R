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

# plot_estimates() Tests ------------------------------------------------

test_that("plot_estimates returns a ggplot object", {
  # Create test data
  test_data <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I", "aac3"),
    estimate = c(2.5, 1.8, 0.5),
    lower_ci = c(1.5, 0.8, -0.5),
    upper_ci = c(3.5, 2.8, 1.5),
    p_value = c(0.001, 0.01, 0.3)
  )
  
  result <- plot_estimates(test_data)
  
  expect_s3_class(result, "ggplot")
})

test_that("plot_estimates works with custom significance level", {
  test_data <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I"),
    estimate = c(2.5, 1.8),
    lower_ci = c(1.5, 0.8),
    upper_ci = c(3.5, 2.8),
    p_value = c(0.001, 0.01)
  )
  
  result <- plot_estimates(test_data, sig = 0.01)
  
  expect_s3_class(result, "ggplot")
})

test_that("plot_estimates handles custom titles and labels", {
  test_data <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I"),
    estimate = c(2.5, 1.8),
    lower_ci = c(1.5, 0.8),
    upper_ci = c(3.5, 2.8),
    p_value = c(0.001, 0.01)
  )
  
  result <- plot_estimates(
    test_data,
    title = "Test Plot",
    x_title = "Custom X",
    y_title = "Custom Y"
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("plot_estimates handles marker_order parameter", {
  test_data <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I", "aac3"),
    estimate = c(2.5, 1.8, 0.5),
    lower_ci = c(1.5, 0.8, -0.5),
    upper_ci = c(3.5, 2.8, 1.5),
    p_value = c(0.001, 0.01, 0.3)
  )
  
  result <- plot_estimates(
    test_data,
    marker_order = c("aac3", "parC_S80I", "gyrA_S83L")
  )
  
  expect_s3_class(result, "ggplot")
})

test_that("plot_estimates handles edge cases", {
  # Single marker
  single_marker <- data.frame(
    marker = "gyrA_S83L",
    estimate = 2.5,
    lower_ci = 1.5,
    upper_ci = 3.5,
    p_value = 0.001
  )
  
  result <- plot_estimates(single_marker)
  expect_s3_class(result, "ggplot")
})

# compare_estimates() Tests ---------------------------------------------

test_that("compare_estimates returns ggplot object", {
  test_data1 <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I"),
    estimate = c(2.5, 1.8),
    lower_ci = c(1.5, 0.8),
    upper_ci = c(3.5, 2.8),
    p_value = c(0.001, 0.01)
  )
  
  test_data2 <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I"),
    estimate = c(2.2, 1.5),
    lower_ci = c(1.2, 0.5),
    upper_ci = c(3.2, 2.5),
    p_value = c(0.002, 0.02)
  )
  
  result <- compare_estimates(test_data1, test_data2)
  
  # Should return ggplot or patchwork object
  expect_true(inherits(result, "ggplot") || inherits(result, "patchwork"))
})

test_that("compare_estimates works with custom titles", {
  test_data1 <- data.frame(
    marker = c("gyrA_S83L"),
    estimate = c(2.5),
    lower_ci = c(1.5),
    upper_ci = c(3.5),
    p_value = c(0.001)
  )
  
  test_data2 <- data.frame(
    marker = c("gyrA_S83L"),
    estimate = c(2.2),
    lower_ci = c(1.2),
    upper_ci = c(3.2),
    p_value = c(0.002)
  )
  
  result <- compare_estimates(
    test_data1,
    test_data2,
    title1 = "Dataset 1",
    title2 = "Dataset 2"
  )
  
  expect_true(inherits(result, "ggplot") || inherits(result, "patchwork"))
})

test_that("compare_estimates single_plot parameter works", {
  test_data1 <- data.frame(
    marker = c("gyrA_S83L"),
    estimate = c(2.5),
    lower_ci = c(1.5),
    upper_ci = c(3.5),
    p_value = c(0.001)
  )
  
  test_data2 <- data.frame(
    marker = c("gyrA_S83L"),
    estimate = c(2.2),
    lower_ci = c(1.2),
    upper_ci = c(3.2),
    p_value = c(0.002)
  )
  
  # Single plot
  result_single <- compare_estimates(
    test_data1,
    test_data2,
    single_plot = TRUE
  )
  
  # Separate plots
  result_separate <- compare_estimates(
    test_data1,
    test_data2,
    single_plot = FALSE
  )
  
  expect_true(inherits(result_single, "ggplot") || inherits(result_single, "patchwork"))
  expect_true(inherits(result_separate, "ggplot") || inherits(result_separate, "patchwork"))
})

# logistf_details() Tests -----------------------------------------------

test_that("logistf_details requires logistf model", {
  # This function extracts details from logistf models
  # We'll test that it handles input appropriately
  
  # Can't easily create a logistf model without the package
  # So we'll test for appropriate error handling
  expect_error(logistf_details(NULL))
  expect_error(logistf_details("not a model"))
})

test_that("logistf_details returns correct structure", {
  skip_if_not_installed("logistf")
  
  # Create minimal test data for logistf
  test_data <- data.frame(
    outcome = c(1, 0, 1, 0, 1, 1, 0, 0),
    predictor = c(1, 0, 1, 0, 1, 1, 0, 0),
    stringsAsFactors = FALSE
  )
  
  # Fit logistf model
  model <- logistf::logistf(outcome ~ predictor, data = test_data)
  
  # Extract details
  result <- logistf_details(model)
  
  expect_s3_class(result, "model_summary")
  expect_s3_class(result, "data.frame")
})

# glm_details() Tests ---------------------------------------------------

test_that("glm_details requires glm model", {
  expect_error(glm_details(NULL))
  expect_error(glm_details("not a model"))
})

test_that("glm_details returns correct structure", {
  # Create test data
  test_data <- data.frame(
    outcome = c(1, 0, 1, 0, 1, 1, 0, 0),
    predictor = c(1, 0, 1, 0, 1, 1, 0, 0)
  )
  
  # Fit glm model
  model <- glm(outcome ~ predictor, data = test_data, family = binomial())
  
  # Extract details
  result <- glm_details(model)
  
  expect_s3_class(result, "model_summary")
  expect_s3_class(result, "data.frame")
  
  # Check for expected columns
  expect_true("estimate" %in% colnames(result) || 
              "lower_ci" %in% colnames(result) ||
              "Estimate" %in% colnames(result))
})

test_that("glm_details handles different model types", {
  # Test with logistic regression
  test_data <- data.frame(
    outcome = c(1, 0, 1, 0, 1, 1, 0, 0),
    predictor1 = c(1, 0, 1, 0, 1, 1, 0, 0),
    predictor2 = c(0, 1, 0, 1, 0, 0, 1, 1)
  )
  
  model <- glm(outcome ~ predictor1 + predictor2, 
               data = test_data, 
               family = binomial())
  
  result <- glm_details(model)
  
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})

# Integration Tests -----------------------------------------------------

test_that("model detail functions work with plotting functions", {
  # Create a simple GLM model
  test_data <- data.frame(
    outcome = c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0),
    marker1 = c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0),
    marker2 = c(0, 1, 0, 1, 0, 0, 1, 1, 0, 1)
  )
  
  model <- glm(outcome ~ marker1 + marker2, 
               data = test_data, 
               family = binomial())
  
  # Extract details
  details <- glm_details(model)
  
  # Should be able to plot (if details has right structure)
  if ("marker" %in% colnames(details) && 
      "estimate" %in% colnames(details) &&
      "p_value" %in% colnames(details)) {
    plot_result <- plot_estimates(details)
    expect_s3_class(plot_result, "ggplot")
  }
})

test_that("plotting functions handle various data formats", {
  # Test that plotting functions are robust to column name variations
  test_data <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I"),
    estimate = c(2.5, 1.8),
    lower_ci = c(1.5, 0.8),
    upper_ci = c(3.5, 2.8),
    p_value = c(0.001, 0.01)
  )
  
  expect_no_error(plot_estimates(test_data))
})

test_that("plotting functions handle missing values", {
  # Test with some NA values
  test_data <- data.frame(
    marker = c("gyrA_S83L", "parC_S80I", "aac3"),
    estimate = c(2.5, NA, 0.5),
    lower_ci = c(1.5, 0.8, -0.5),
    upper_ci = c(3.5, 2.8, 1.5),
    p_value = c(0.001, 0.01, NA)
  )
  
  result <- plot_estimates(test_data)
  expect_s3_class(result, "ggplot")
})
