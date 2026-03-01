# AMRgen Test Suite

This directory contains comprehensive unit tests for the AMRgen package using testthat v3.

## Test Files Overview

### Data Tests
- **test-data.R** - Validates package data objects (ecoli_ast_raw, ecoli_ast, ecoli_geno_raw)

### Import/Export Tests
- **test-import_pheno.R** - Tests import functions for AST data from multiple sources (NCBI, EBI, VITEK, etc.)
- **test-download_ncbi_ast.R** - Tests NCBI AST data download functionality (requires internet)
- **test-gtdb.R** - Tests GTDB species name parsing and import

### Analysis Tests
- **test-breakpoints.R** - Tests clinical breakpoint retrieval with hierarchical fallback
- **test-get_binary_matrix.R** - Tests genotype-phenotype binary matrix generation
- **test-eucast_distributions.R** - Tests EUCAST wild-type distribution retrieval (requires internet)

### Visualization Tests
- **test-plot_estimates.R** - Tests logistic regression coefficient plotting
- **test-assay_distribution.R** - Tests MIC/disk distribution visualization
- **test-amr_upset.R** - Tests upset plot generation for AMR marker combinations

### Integration Tests
- **test-hAMRonization.R** - Tests data harmonization (requires Python/reticulate)

## Running Tests

### Run all tests
```r
devtools::test()
```

### Run specific test file
```r
testthat::test_file("tests/testthat/test-data.R")
```

### Run tests with coverage
```r
covr::package_coverage()
```

## Test Dependencies

Some tests require:
- **Internet connection**: Tests with `skip_if_offline()`
- **Python/reticulate**: hAMRonization tests with `skip_if_not_installed("reticulate")`
- **External packages**: logistf for some model tests

Tests are designed to skip gracefully when dependencies are unavailable.

## Test Data

Tests use package example datasets:
- `ecoli_ast` - Processed E. coli AST data
- `ecoli_ast_raw` - Raw NCBI AST data
- `ecoli_geno_raw` - AMRFinderPlus genotype data

## CI/CD

Tests run automatically via GitHub Actions on:
- Multiple R versions (devel, release, oldrel)
- Multiple platforms (macOS, Ubuntu, Windows)
- Pull requests and pushes to any branch
