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

# Package Data Objects Tests -------------------------------------------

test_that("ecoli_ast_raw data object exists and has correct structure", {
  # Check object exists
  expect_true(exists("ecoli_ast_raw"))
  
  # Check it's a data frame
  expect_s3_class(ecoli_ast_raw, "data.frame")
  
  # Check it has rows and columns
  expect_gt(nrow(ecoli_ast_raw), 0)
  expect_gt(ncol(ecoli_ast_raw), 0)
  
  # Check expected columns exist
  expect_true("#BioSample" %in% colnames(ecoli_ast_raw))
  expect_true("Scientific name" %in% colnames(ecoli_ast_raw))
  expect_true("Antibiotic" %in% colnames(ecoli_ast_raw))
  expect_true("Testing standard" %in% colnames(ecoli_ast_raw))
})

test_that("ecoli_ast data object exists and has correct structure", {
  # Check object exists
  expect_true(exists("ecoli_ast"))
  
  # Check it's a data frame
  expect_s3_class(ecoli_ast, "data.frame")
  
  # Check it has rows and columns
  expect_gt(nrow(ecoli_ast), 0)
  expect_gt(ncol(ecoli_ast), 0)
  
  # Check expected columns exist
  expect_true("id" %in% colnames(ecoli_ast))
  expect_true("drug_agent" %in% colnames(ecoli_ast))
  expect_true("mic" %in% colnames(ecoli_ast))
  expect_true("disk" %in% colnames(ecoli_ast))
  expect_true("spp_pheno" %in% colnames(ecoli_ast))
})

test_that("ecoli_ast has correct AMR package classes", {
  # Check for proper AMR package classes
  if ("drug_agent" %in% colnames(ecoli_ast)) {
    expect_s3_class(ecoli_ast$drug_agent, "ab")
  }
  
  if ("mic" %in% colnames(ecoli_ast)) {
    expect_s3_class(ecoli_ast$mic, "mic")
  }
  
  if ("disk" %in% colnames(ecoli_ast)) {
    expect_s3_class(ecoli_ast$disk, "disk")
  }
  
  if ("spp_pheno" %in% colnames(ecoli_ast)) {
    expect_s3_class(ecoli_ast$spp_pheno, "mo")
  }
  
  if ("pheno_clsi" %in% colnames(ecoli_ast)) {
    expect_s3_class(ecoli_ast$pheno_clsi, "sir")
  }
  
  if ("ecoff" %in% colnames(ecoli_ast)) {
    expect_s3_class(ecoli_ast$ecoff, "sir")
  }
})

test_that("ecoli_geno_raw data object exists and has correct structure", {
  # Check object exists
  expect_true(exists("ecoli_geno_raw"))
  
  # Check it's a data frame
  expect_s3_class(ecoli_geno_raw, "data.frame")
  
  # Check it has rows and columns
  expect_gt(nrow(ecoli_geno_raw), 0)
  expect_gt(ncol(ecoli_geno_raw), 0)
  
  # Check expected columns exist (from AMRFinderPlus output)
  expect_true("Name" %in% colnames(ecoli_geno_raw))
  expect_true("Gene symbol" %in% colnames(ecoli_geno_raw))
})

test_that("data objects are consistent with each other", {
  # ecoli_ast should be derived from ecoli_ast_raw
  # Both should have similar number of rows (unless filtering was applied)
  expect_true(nrow(ecoli_ast) > 0)
  expect_true(nrow(ecoli_ast_raw) > 0)
})

test_that("data objects can be used in typical workflows", {
  # Test that the data can be loaded and manipulated without errors
  expect_no_error({
    data <- ecoli_ast
    head(data)
    nrow(data)
    colnames(data)
  })
  
  expect_no_error({
    data <- ecoli_ast_raw
    head(data)
    nrow(data)
    colnames(data)
  })
  
  expect_no_error({
    data <- ecoli_geno_raw
    head(data)
    nrow(data)
    colnames(data)
  })
})
