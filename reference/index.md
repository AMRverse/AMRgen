# Package index

## Data Import & Preprocessing

Functions to import, harmonise, and prepare genotype and phenotype data
from public repositories or internal formats.

- [`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
  : Import and Process AMRFinderPlus Results
- [`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md)
  : Import and Process AST Data from an EBI or NCBI antibiogram File
- [`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md)
  : Import and Process AST Data from files downloaded from the EBI AMR
  portal website
- [`import_ebi_ast_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md)
  : Import and Process AST Data files retrieved from the EBI AMR portal
  FTP site
- [`import_gtdb()`](https://AMRverse.github.io/AMRgen/reference/import_gtdb.md)
  : Import GTDB Output
- [`import_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md)
  : Import and Process AST Data from an NCBI File
- [`import_vitek_ast()`](https://AMRverse.github.io/AMRgen/reference/import_vitek_ast.md)
  : Import and Process AST Data from Vitek Output Files
- [`import_whonet_ast()`](https://AMRverse.github.io/AMRgen/reference/import_whonet_ast.md)
  : Import and Process AST Data from WHONET Output Files
- [`export_ast()`](https://AMRverse.github.io/AMRgen/reference/export_ast.md)
  : Export AST Data
- [`export_ebi_antibiogram()`](https://AMRverse.github.io/AMRgen/reference/export_ebi_antibiogram.md)
  : Export EBI Antibiogram
- [`export_ncbi_biosample()`](https://AMRverse.github.io/AMRgen/reference/export_ncbi_biosample.md)
  : Export NCBI BioSample Antibiogram
- [`convert_aa_code()`](https://AMRverse.github.io/AMRgen/reference/convert_aa_code.md)
  : Convert single-letter amino acid code(s) to three-letter code(s)
- [`convert_mutation()`](https://AMRverse.github.io/AMRgen/reference/convert_mutation.md)
  : Convert mutation string based on method
- [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md)
  : Get Binary Matrix of Genotype and Phenotype Data
- [`get_combo_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_combo_matrix.md)
  : Generate matrix of marker combinations
- [`as.gene()`](https://AMRverse.github.io/AMRgen/reference/as.gene.md)
  : Gene Class and AMR Parsing Functions
- [`gtdb.mo()`](https://AMRverse.github.io/AMRgen/reference/gtdb.mo.md)
  : Get Microorganism from GTDB Species Name
- [`format_ast()`](https://AMRverse.github.io/AMRgen/reference/format_ast.md)
  : Import and Process AST Data from a generic format
- [`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md)
  : Download antimicrobial genotype or phenotype data from the EBI AMR
  Portal

## Resistance Interpretation

Core tools for interpreting antimicrobial resistance based on EUCAST
breakpoints and custom models.

- [`compare_estimates()`](https://AMRverse.github.io/AMRgen/reference/compare_estimates.md)
  : Plot to Compare Two Sets of Estimates
- [`compare_geno_pheno_id()`](https://AMRverse.github.io/AMRgen/reference/compare_geno_pheno_id.md)
  : Compare Genotype and Phenotype Data by Sample ID
- [`get_eucast_amr_distribution()`](https://AMRverse.github.io/AMRgen/reference/get_eucast_amr_distribution.md)
  [`get_eucast_mic_distribution()`](https://AMRverse.github.io/AMRgen/reference/get_eucast_amr_distribution.md)
  [`get_eucast_disk_distribution()`](https://AMRverse.github.io/AMRgen/reference/get_eucast_amr_distribution.md)
  [`compare_mic_with_eucast()`](https://AMRverse.github.io/AMRgen/reference/get_eucast_amr_distribution.md)
  [`compare_disk_with_eucast()`](https://AMRverse.github.io/AMRgen/reference/get_eucast_amr_distribution.md)
  : Get and Compare Antimicrobial Wild Type Distributions from EUCAST
- [`eucast_supported_ab_distributions()`](https://AMRverse.github.io/AMRgen/reference/eucast_supported_ab_distributions.md)
  : Retrieve Available Antimicrobial Wild Type Distributions from EUCAST
- [`merge_logreg_soloppv()`](https://AMRverse.github.io/AMRgen/reference/merge_logreg_soloppv.md)
  : Merge Logistic Regression and Solo PPV Statistics
- [`interpret_ast()`](https://AMRverse.github.io/AMRgen/reference/interpret_ast.md)
  : Interpret AST data in a standard format tibble

## Modelling and analysis

Statistical models for resistance prediction and inference, including
logistic regression and Firth regression.

- [`solo_ppv_analysis()`](https://AMRverse.github.io/AMRgen/reference/solo_ppv_analysis.md)
  : Perform Solo PPV Analysis for AMR Markers
- [`amr_logistic()`](https://AMRverse.github.io/AMRgen/reference/amr_logistic.md)
  : AMR Logistic Regression Analysis
- [`glm_details()`](https://AMRverse.github.io/AMRgen/reference/glm_details.md)
  : Extract Details from a Generalized Linear Model
- [`logistf_details()`](https://AMRverse.github.io/AMRgen/reference/logistf_details.md)
  : Extract Details from a logistf Model
- [`getBreakpoints()`](https://AMRverse.github.io/AMRgen/reference/getBreakpoints.md)
  : Get Clinical Breakpoints for an Antibiotic
- [`checkBreakpoints()`](https://AMRverse.github.io/AMRgen/reference/checkBreakpoints.md)
  : Check and Retrieve Breakpoints for an Antibiotic

## Visualisation & Reporting

Tools to visualise results, including MIC distributions, model outputs,
and genotype-phenotype relationships.

- [`plot_combined_stats()`](https://AMRverse.github.io/AMRgen/reference/plot_combined_stats.md)
  : Plot Combined Statistics
- [`plot_estimates()`](https://AMRverse.github.io/AMRgen/reference/plot_estimates.md)
  : Plot Estimates from a Table of Results
- [`plot_solo_logReg()`](https://AMRverse.github.io/AMRgen/reference/plot_solo_logReg.md)
  : Plot Combined Statistics of Logistic Regression and Solo PPV
- [`amr_upset()`](https://AMRverse.github.io/AMRgen/reference/amr_upset.md)
  : Generate a Series of Plots for AMR Gene and Combination Analysis
- [`assay_by_var()`](https://AMRverse.github.io/AMRgen/reference/assay_by_var.md)
  : Generate a Stacked Bar Plot of Assay Values Colored by a Variable

## Data

Example datasets for demonstration and reproducible analysis.

- [`ecoli_ast`](https://AMRverse.github.io/AMRgen/reference/ecoli_ast.md)
  : E. coli NCBI AST Example Data, Re-interpreted with AMR Package
- [`ecoli_ast_raw`](https://AMRverse.github.io/AMRgen/reference/ecoli_ast_raw.md)
  : E. coli NCBI AST Example Data
- [`ecoli_geno_raw`](https://AMRverse.github.io/AMRgen/reference/ecoli_geno_raw.md)
  : E. coli Genotype Example Data
