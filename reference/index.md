# Package index

## Data Import & Preprocessing

Functions to import, harmonise, and prepare genotype and phenotype data
from public repositories or internal formats.

- [`import_abricate()`](https://AMRverse.github.io/AMRgen/reference/import_abricate.md)
  : Import and Process Abricate Results
- [`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
  : Import and Process AMRFinderPlus Results
- [`import_amrfp_ebi()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp_ebi.md)
  : Import EBI-processed AMRFinderPlus Genotypes
- [`import_amrfp_ebi_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp_ebi_ftp.md)
  : Import EBI-processed AMRFinderPlus Genotypes from FTP
- [`import_amrfp_ebi_web()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp_ebi_web.md)
  : Import EBI-processed AMRFinderPlus Genotypes from Web
- [`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md)
  : Import and process antimicrobial phenotype data from common sources
- [`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md)
  : Import and process antimicrobial susceptibility phenotype data from
  the EBI AMR web portal
- [`import_ebi_ast_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md)
  : Import and process antimicrobial phenotype data files retrieved from
  the EBI AMR portal FTP site
- [`import_geno()`](https://AMRverse.github.io/AMRgen/reference/import_geno.md)
  : Import and process antimicrobial genotype data from common sources
- [`import_gtdb()`](https://AMRverse.github.io/AMRgen/reference/import_gtdb.md)
  : Import GTDB Output
- [`import_kleborate()`](https://AMRverse.github.io/AMRgen/reference/import_kleborate.md)
  : Import and Process Kleborate Results
- [`import_microscan_ast()`](https://AMRverse.github.io/AMRgen/reference/import_microscan_ast.md)
  : Import and process antimicrobial phenotype data exported from
  MicroScan instruments
- [`import_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md)
  : Import and process antimicrobial susceptibility phenotype data from
  the NCBI AST browser
- [`import_ncbi_biosample()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_biosample.md)
  : Import and process antimicrobial phenotype data retrieved from NCBI
  BioSamples
- [`import_pheno()`](https://AMRverse.github.io/AMRgen/reference/import_pheno.md)
  : Import and process antimicrobial phenotype data from common sources
- [`import_phoenix_ast()`](https://AMRverse.github.io/AMRgen/reference/import_phoenix_ast.md)
  : Import and process antimicrobial phenotype data exported from BD
  Phoenix instruments
- [`import_sensititre_ast()`](https://AMRverse.github.io/AMRgen/reference/import_sensititre_ast.md)
  : Import and process antimicrobial phenotype data exported from
  Sensititre instruments
- [`import_vitek_ast()`](https://AMRverse.github.io/AMRgen/reference/import_vitek_ast.md)
  : Import and process antimicrobial phenotype data exported from Vitek
  instruments
- [`import_whonet_ast()`](https://AMRverse.github.io/AMRgen/reference/import_whonet_ast.md)
  : Import and process antimicrobial phenotype data from WHONET files
- [`export_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/export_ebi_ast.md)
  : Export EBI Antibiogram
- [`export_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/export_ncbi_ast.md)
  : Export NCBI BioSample Antibiogram
- [`convert_aa_code()`](https://AMRverse.github.io/AMRgen/reference/convert_aa_code.md)
  : Convert single-letter amino acid code(s) to three-letter code(s)
- [`convert_mutation()`](https://AMRverse.github.io/AMRgen/reference/convert_mutation.md)
  : Convert mutation string based on method
- [`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md)
  : Download antimicrobial genotype or phenotype data from the EBI AMR
  Portal
- [`download_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/download_ncbi_ast.md)
  : Download NCBI antimicrobial susceptibility testing (AST) data
- [`summarise_geno()`](https://AMRverse.github.io/AMRgen/reference/summarise_geno.md)
  : Summarise a Genotype Table
- [`summarise_pheno()`](https://AMRverse.github.io/AMRgen/reference/summarise_pheno.md)
  : Summarise a Phenotype Table
- [`format_ast()`](https://AMRverse.github.io/AMRgen/reference/format_ast.md)
  : Import and process antimicrobial phenotype data from a generic
  format
- [`format_ebi_json()`](https://AMRverse.github.io/AMRgen/reference/format_ebi_json.md)
  : Generate EBI antibiogram submission in JSON
- [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md)
  : Get Binary Matrix of Genotype and Phenotype Data
- [`get_combo_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_combo_matrix.md)
  : Add marker combinations to a binary geno-pheno matrix
- [`combo_stats()`](https://AMRverse.github.io/AMRgen/reference/combo_stats.md)
  : Generate a Series of Plots for AMR Gene and Combination Analysis
- [`as.gene()`](https://AMRverse.github.io/AMRgen/reference/as.gene.md)
  : Gene Class and AMR Parsing Functions
- [`gtdb.mo()`](https://AMRverse.github.io/AMRgen/reference/gtdb.mo.md)
  : Get Microorganism from GTDB Species Name
- [`query_ncbi_bq_ast()`](https://AMRverse.github.io/AMRgen/reference/query_ncbi_bq_ast.md)
  : Query antimicrobial phenotype (AST) data from NCBI Pathogen
  Detection BigQuery
- [`query_ncbi_bq_geno()`](https://AMRverse.github.io/AMRgen/reference/query_ncbi_bq_geno.md)
  : Query antimicrobial genotype (MicroBIGG-E) data from NCBI Pathogen
  Detection BigQuery
- [`hamronize_data()`](https://AMRverse.github.io/AMRgen/reference/hamronize_data.md)
  : hamronize_data

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
  : Interpret antimicrobial susceptibility phenotype data in a standard
  format tibble

## Modelling and analysis

Statistical models for resistance prediction and inference, including
logistic regression and Firth regression.

- [`amr_upset()`](https://AMRverse.github.io/AMRgen/reference/amr_upset.md)
  : Generate Upset Plot
- [`solo_ppv_analysis()`](https://AMRverse.github.io/AMRgen/reference/solo_ppv_analysis.md)
  : Perform Solo PPV Analysis for AMR Markers
- [`ppv()`](https://AMRverse.github.io/AMRgen/reference/ppv.md) :
  Generate PPV Plot
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
- [`concordance()`](https://AMRverse.github.io/AMRgen/reference/concordance.md)
  : Calculate genotype-phenotype concordance from binary matrix
- [`print(`*`<amr_concordance>`*`)`](https://AMRverse.github.io/AMRgen/reference/print.amr_concordance.md)
  : Print method for amr_concordance objects

## Visualisation & Reporting

Tools to visualise results, including MIC distributions, model outputs,
and genotype-phenotype relationships.

- [`plot_combined_stats()`](https://AMRverse.github.io/AMRgen/reference/plot_combined_stats.md)
  : Plot Combined Statistics
- [`plot_estimates()`](https://AMRverse.github.io/AMRgen/reference/plot_estimates.md)
  : Plot Estimates from a Table of Results
- [`plot_solo_logReg()`](https://AMRverse.github.io/AMRgen/reference/plot_solo_logReg.md)
  : Plot Combined Statistics of Logistic Regression and Solo PPV
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

- [`ecoli_cip_vs_ref`](https://AMRverse.github.io/AMRgen/reference/ecoli_cip_vs_ref.md)
  : E. coli Ciprofloxacin MIC Distribution Example Data

- [`ecoli_cip_mic_data`](https://AMRverse.github.io/AMRgen/reference/ecoli_cip_mic_data.md)
  : EUCAST Reference distribution for Ciprofloxacin in E. coli

- [`staph_ast_ncbi`](https://AMRverse.github.io/AMRgen/reference/staph_ast_ncbi.md)
  : S. aureus Example of Imported NCBI Phenotype Data

- [`staph_ast_ncbi_raw`](https://AMRverse.github.io/AMRgen/reference/staph_ast_ncbi_raw.md)
  : S. aureus Example of Raw Phenotype Data Downloaded from NCBI
  BioSamples via Entrez API

- [`staph_ast_ncbi_cloud_raw`](https://AMRverse.github.io/AMRgen/reference/staph_ast_ncbi_cloud_raw.md)
  : S. aureus Example of Raw Phenotype Data Downloaded from NCBI via
  Google Cloud BigQuery

- [`staph_geno_ncbi_cloud_raw`](https://AMRverse.github.io/AMRgen/reference/staph_geno_ncbi_cloud_raw.md)
  : S. aureus Example of Raw Genotype Data Downloaded from NCBI via
  Google Cloud BigQuery

- [`staph_ast_ebi`](https://AMRverse.github.io/AMRgen/reference/staph_ast_ebi.md)
  : S. aureus Example of Imported EBI Phenotype Data

- [`staph_geno_ebi`](https://AMRverse.github.io/AMRgen/reference/staph_geno_ebi.md)
  : S. aureus Example of Imported EBI Genotype Data

- [`amrfp_drugs_table`](https://AMRverse.github.io/AMRgen/reference/amrfp_drugs_table.md)
  : NCBI Subclass mapping to drug class

- [`pheno_eco_2075`](https://AMRverse.github.io/AMRgen/reference/pheno_eco_2075.md)
  : E. coli AST data from Mills et al 2022

- [`geno_eco_2075`](https://AMRverse.github.io/AMRgen/reference/geno_eco_2075.md)
  : E. coli genotype data from Mills et al 2022

- [`eurogasp_geno_raw`](https://AMRverse.github.io/AMRgen/reference/eurogasp_geno_raw.md)
  : N. gonorrhoeae Euro-GASP Genotype Data Use Case 1

- [`ngono_cro_geno_raw`](https://AMRverse.github.io/AMRgen/reference/ngono_cro_geno_raw.md)
  : N. gonorrhoeae PBP2 Mutations Genotype Data Use Case 2

- [`ngono_tet_geno_raw`](https://AMRverse.github.io/AMRgen/reference/ngono_tet_geno_raw.md)
  : N. gonorrhoeae Tetracycline Resistance Genotype Data Use Case 3

- [`ngono_cro_pheno_raw`](https://AMRverse.github.io/AMRgen/reference/ngono_cro_pheno_raw.md)
  : N. gonorrhoeae PBP2 Mutations Phenotype Data Use Case 2

- [`eurogasp_pheno_raw`](https://AMRverse.github.io/AMRgen/reference/eurogasp_pheno_raw.md)
  : N. gonorrhoeae Euro-GASP Phenotype Data Use Case 1

- [`ngono_tet_pheno_raw`](https://AMRverse.github.io/AMRgen/reference/ngono_tet_pheno_raw.md)
  : N. gonorrhoeae Tetracycline Resistance Phenotype Data Use Case 3

- [`salm_raw`](https://AMRverse.github.io/AMRgen/reference/salm_raw.md)
  : Example Salmonella Genotype-Phenotype Data

- [`kleborate_raw`](https://AMRverse.github.io/AMRgen/reference/kleborate_raw.md)
  : Example Kleborate Genotype Data from EuSCAPE project

- [`kleborate_raw_v313`](https://AMRverse.github.io/AMRgen/reference/kleborate_raw_v313.md)
  : Example Kleborate v3.1.3 Genotype Data from EuSCAPE project

- [`kleborate_classes`](https://AMRverse.github.io/AMRgen/reference/kleborate_classes.md)
  : Table mapping Kleborate drug class columns

- [`ast_CLI_public`](https://AMRverse.github.io/AMRgen/reference/ast_CLI_public.md)
  :

  Clindamycin MIC data for 5914 *Staphylococcus aureus* isolates

- [`afp_CLI_public`](https://AMRverse.github.io/AMRgen/reference/afp_CLI_public.md)
  : S. aureus Clindamycin Resistance Genotype Data

- [`ST_data_CLI`](https://AMRverse.github.io/AMRgen/reference/ST_data_CLI.md)
  :

  ST data for *Staphylococcus aureus* genomes for clindamycin vignette
