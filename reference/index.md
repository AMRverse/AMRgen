# Package index

## Data Import & Preprocessing

Functions to import, harmonise, and prepare genotype and phenotype data
from public repositories or internal formats.

- [`import_abricate()`](https://amrgen.org/reference/import_abricate.md)
  : Import and Process ABRicate Results
- [`import_amrfp()`](https://amrgen.org/reference/import_amrfp.md) :
  Import and Process AMRFinderPlus Results
- [`import_amrfp_ebi()`](https://amrgen.org/reference/import_amrfp_ebi.md)
  : Import EBI-processed AMRFinderPlus Genotypes
- [`import_amrfp_ebi_ftp()`](https://amrgen.org/reference/import_amrfp_ebi_ftp.md)
  : Import EBI-processed AMRFinderPlus Genotypes from FTP
- [`import_amrfp_ebi_web()`](https://amrgen.org/reference/import_amrfp_ebi_web.md)
  : Import EBI-processed AMRFinderPlus Genotypes from Web
- [`import_amrrules_predictions()`](https://amrgen.org/reference/import_amrrules_predictions.md)
  : Import phenotype predictions from AMRrules output
- [`import_ebi_pheno()`](https://amrgen.org/reference/import_ebi_pheno.md)
  : Import and process antimicrobial susceptibility phenotype data from
  the EBI AMR web portal
- [`import_ebi_pheno_ftp()`](https://amrgen.org/reference/import_ebi_pheno_ftp.md)
  : Import and process antimicrobial phenotype data files retrieved from
  the EBI AMR portal FTP site
- [`import_geno()`](https://amrgen.org/reference/import_geno.md) :
  Import and process antimicrobial genotype data from common sources
- [`import_kleborate()`](https://amrgen.org/reference/import_kleborate.md)
  : Import and Process Kleborate Results
- [`import_microscan_pheno()`](https://amrgen.org/reference/import_microscan_pheno.md)
  : Import and process antimicrobial phenotype data exported from
  MicroScan instruments
- [`import_ncbi_biosample()`](https://amrgen.org/reference/import_ncbi_biosample.md)
  : Import and process antimicrobial phenotype data retrieved from NCBI
  BioSamples
- [`import_ncbi_pheno()`](https://amrgen.org/reference/import_ncbi_pheno.md)
  : Import and process antimicrobial susceptibility phenotype data from
  the NCBI AST browser
- [`import_pheno()`](https://amrgen.org/reference/import_pheno.md) :
  Import and process antimicrobial phenotype data from common sources
- [`import_phoenix_pheno()`](https://amrgen.org/reference/import_phoenix_pheno.md)
  : Import and process antimicrobial phenotype data exported from BD
  Phoenix instruments
- [`import_rgi()`](https://amrgen.org/reference/import_rgi.md) : Import
  and Process Resistance Gene Identifier (RGI) Results
- [`import_sensititre_pheno()`](https://amrgen.org/reference/import_sensititre_pheno.md)
  : Import and process antimicrobial phenotype data exported from
  Sensititre instruments
- [`import_sirscan_pheno()`](https://amrgen.org/reference/import_sirscan_pheno.md)
  : Import and process antimicrobial susceptibility data from SIRscan
  (Bio-Rad)
- [`import_vitek_pheno()`](https://amrgen.org/reference/import_vitek_pheno.md)
  : Import and process antimicrobial phenotype data exported from Vitek
  instruments
- [`import_whonet_pheno()`](https://amrgen.org/reference/import_whonet_pheno.md)
  : Import and process antimicrobial phenotype data from WHONET files
- [`export_ebi_pheno()`](https://amrgen.org/reference/export_ebi_pheno.md)
  : Export EBI Antibiogram
- [`export_ncbi_pheno()`](https://amrgen.org/reference/export_ncbi_pheno.md)
  : Export NCBI BioSample Antibiogram
- [`convert_aa_code()`](https://amrgen.org/reference/convert_aa_code.md)
  : Convert single-letter amino acid code(s) to three-letter code(s)
- [`convert_mutation()`](https://amrgen.org/reference/convert_mutation.md)
  : Convert mutation string based on method
- [`download_ebi()`](https://amrgen.org/reference/download_ebi.md) :
  Download antimicrobial genotype or phenotype data from the EBI AMR
  Portal
- [`download_ncbi_pheno()`](https://amrgen.org/reference/download_ncbi_pheno.md)
  : Download NCBI antimicrobial susceptibility testing (AST) data
- [`summarise_geno()`](https://amrgen.org/reference/summarise_geno.md) :
  Summarise a Genotype Table
- [`summarise_geno_pheno()`](https://amrgen.org/reference/summarise_geno_pheno.md)
  : Summarise the intersection of a genotype table and a phenotype table
- [`summarise_pheno()`](https://amrgen.org/reference/summarise_pheno.md)
  : Summarise a Phenotype Table
- [`format_ebi_json()`](https://amrgen.org/reference/format_ebi_json.md)
  : Generate EBI antibiogram submission in JSON
- [`format_pheno()`](https://amrgen.org/reference/format_pheno.md) :
  Import and process antimicrobial phenotype data from a generic format
- [`get_binary_matrix()`](https://amrgen.org/reference/get_binary_matrix.md)
  : Get Binary Matrix of Genotype and Phenotype Data
- [`get_combo_matrix()`](https://amrgen.org/reference/get_combo_matrix.md)
  : Add marker combinations to a binary geno-pheno matrix
- [`query_ncbi_bq_geno()`](https://amrgen.org/reference/query_ncbi_bq_geno.md)
  : Query antimicrobial genotype (MicroBIGG-E) data from NCBI Pathogen
  Detection BigQuery
- [`query_ncbi_bq_pheno()`](https://amrgen.org/reference/query_ncbi_bq_pheno.md)
  : Query antimicrobial phenotype data from NCBI Pathogen Detection
  BigQuery

## Resistance Interpretation

Core tools for interpreting antimicrobial resistance based on EUCAST
breakpoints and custom models.

- [`compare_estimates()`](https://amrgen.org/reference/compare_estimates.md)
  : Plot to Compare Two Sets of Estimates
- [`compare_geno_pheno_id()`](https://amrgen.org/reference/compare_geno_pheno_id.md)
  : Compare Genotype and Phenotype Data by Sample ID
- [`get_eucast_amr_distribution()`](https://amrgen.org/reference/get_eucast_amr_distribution.md)
  [`get_eucast_mic_distribution()`](https://amrgen.org/reference/get_eucast_amr_distribution.md)
  [`get_eucast_disk_distribution()`](https://amrgen.org/reference/get_eucast_amr_distribution.md)
  [`compare_mic_with_eucast()`](https://amrgen.org/reference/get_eucast_amr_distribution.md)
  [`compare_disk_with_eucast()`](https://amrgen.org/reference/get_eucast_amr_distribution.md)
  : Get and Compare Antimicrobial Wild Type Distributions from EUCAST
- [`eucast_supported_ab_distributions()`](https://amrgen.org/reference/eucast_supported_ab_distributions.md)
  : Retrieve Available Antimicrobial Wild Type Distributions from EUCAST
- [`merge_logreg_soloppv()`](https://amrgen.org/reference/merge_logreg_soloppv.md)
  : Merge Logistic Regression and Solo PPV Statistics
- [`interpret_pheno()`](https://amrgen.org/reference/interpret_pheno.md)
  : Interpret antimicrobial susceptibility phenotype data in a standard
  format tibble

## Modelling and analysis

Statistical models for resistance prediction and inference, including
logistic regression and Firth regression.

- [`amr_upset()`](https://amrgen.org/reference/amr_upset.md) : Generate
  Upset Plot
- [`solo_ppv()`](https://amrgen.org/reference/solo_ppv.md) : Perform
  Solo PPV Analysis for AMR Markers
- [`ppv()`](https://amrgen.org/reference/ppv.md) : Generate PPV Plot
- [`amr_logistic()`](https://amrgen.org/reference/amr_logistic.md) : AMR
  Logistic Regression Analysis
- [`glm_details()`](https://amrgen.org/reference/glm_details.md) :
  Extract Details from a Generalised Linear Model
- [`logistf_details()`](https://amrgen.org/reference/logistf_details.md)
  : Extract Details from a logistf Model
- [`getBreakpoints()`](https://amrgen.org/reference/getBreakpoints.md) :
  Get Clinical Breakpoints for an Antibiotic
- [`checkBreakpoints()`](https://amrgen.org/reference/checkBreakpoints.md)
  : Check and Retrieve Breakpoints for an Antibiotic
- [`concordance()`](https://amrgen.org/reference/concordance.md) :
  Calculate genotype-phenotype concordance from binary matrix
- [`concordance_from_tables()`](https://amrgen.org/reference/concordance_from_tables.md)
  : Assess concordance between tables of observed and predicted
  phenotypes
- [`print(`*`<amr_concordance>`*`)`](https://amrgen.org/reference/print.amr_concordance.md)
  : Print method for amr_concordance objects

## Visualisation & Reporting

Tools to visualise results, including MIC distributions, model outputs,
and genotype-phenotype relationships.

- [`plot_combined_stats()`](https://amrgen.org/reference/plot_combined_stats.md)
  : Plot Combined Statistics
- [`plot_estimates()`](https://amrgen.org/reference/plot_estimates.md) :
  Plot Estimates from a Table of Results
- [`plot_solo_logReg()`](https://amrgen.org/reference/plot_solo_logReg.md)
  : Plot Combined Statistics of Logistic Regression and Solo PPV
- [`assay_by_var()`](https://amrgen.org/reference/assay_by_var.md) :
  Generate a Stacked Bar Plot of Assay Values Colored by a Variable

## Data

Example datasets for demonstration and reproducible analysis.

- [`ecoli_pheno`](https://amrgen.org/reference/ecoli_pheno.md) : E. coli
  NCBI AST Example Data, Re-interpreted with AMR Package

- [`ecoli_pheno_raw`](https://amrgen.org/reference/ecoli_pheno_raw.md)
  : E. coli NCBI AST Example Data

- [`ecoli_geno_raw`](https://amrgen.org/reference/ecoli_geno_raw.md)
  : E. coli Genotype Example Data

- [`ecoli_cip_vs_ref`](https://amrgen.org/reference/ecoli_cip_vs_ref.md)
  : E. coli Ciprofloxacin MIC Distribution Example Data

- [`ecoli_cip_mic_data`](https://amrgen.org/reference/ecoli_cip_mic_data.md)
  : EUCAST Reference distribution for Ciprofloxacin in E. coli

- [`staph_pheno_ncbi`](https://amrgen.org/reference/staph_pheno_ncbi.md)
  : S. aureus Example of Imported NCBI Phenotype Data

- [`staph_pheno_ncbi_raw`](https://amrgen.org/reference/staph_pheno_ncbi_raw.md)
  : S. aureus Example of Raw Phenotype Data Downloaded from NCBI
  BioSamples via Entrez API

- [`staph_pheno_ncbi_cloud_raw`](https://amrgen.org/reference/staph_pheno_ncbi_cloud_raw.md)
  : S. aureus Example of Raw Phenotype Data Downloaded from NCBI via
  Google Cloud BigQuery

- [`staph_geno_ncbi_cloud_raw`](https://amrgen.org/reference/staph_geno_ncbi_cloud_raw.md)
  : S. aureus Example of Raw Genotype Data Downloaded from NCBI via
  Google Cloud BigQuery

- [`staph_pheno_ebi`](https://amrgen.org/reference/staph_pheno_ebi.md)
  : S. aureus Example of Imported EBI Phenotype Data

- [`staph_geno_ebi`](https://amrgen.org/reference/staph_geno_ebi.md)
  : S. aureus Example of Imported EBI Genotype Data

- [`amrfp_drugs_table`](https://amrgen.org/reference/amrfp_drugs_table.md)
  : NCBI Subclass mapping to drug class

- [`pheno_eco_2075`](https://amrgen.org/reference/pheno_eco_2075.md)
  : E. coli AST data from Mills et al 2022

- [`geno_eco_2075`](https://amrgen.org/reference/geno_eco_2075.md) : E.
  coli genotype data from Mills et al 2022

- [`eurogasp_geno_raw`](https://amrgen.org/reference/eurogasp_geno_raw.md)
  : N. gonorrhoeae Euro-GASP Genotype Data Use Case 1

- [`ngono_cro_geno_raw`](https://amrgen.org/reference/ngono_cro_geno_raw.md)
  : N. gonorrhoeae PBP2 Mutations Genotype Data Use Case 2

- [`ngono_tet_geno_raw`](https://amrgen.org/reference/ngono_tet_geno_raw.md)
  : N. gonorrhoeae Tetracycline Resistance Genotype Data Use Case 3

- [`ngono_cro_pheno_raw`](https://amrgen.org/reference/ngono_cro_pheno_raw.md)
  : N. gonorrhoeae PBP2 Mutations Phenotype Data Use Case 2

- [`eurogasp_pheno_raw`](https://amrgen.org/reference/eurogasp_pheno_raw.md)
  : N. gonorrhoeae Euro-GASP Phenotype Data Use Case 1

- [`azm_comparison`](https://amrgen.org/reference/azm_comparison.md)
  : N. gonorrhoeae Euro-GASP azithromycin data vs EUCAST reference
  distribution

- [`cip_comparison`](https://amrgen.org/reference/cip_comparison.md)
  : N. gonorrhoeae Euro-GASP ciprofloxacin data vs EUCAST reference
  distribution

- [`cro_comparison`](https://amrgen.org/reference/cro_comparison.md)
  : N. gonorrhoeae Euro-GASP ceftriaxone data vs EUCAST reference
  distribution

- [`cfm_comparison`](https://amrgen.org/reference/cfm_comparison.md)
  : N. gonorrhoeae Euro-GASP cefixime data vs EUCAST reference
  distribution

- [`ngono_tet_pheno_raw`](https://amrgen.org/reference/ngono_tet_pheno_raw.md)
  : N. gonorrhoeae Tetracycline Resistance Phenotype Data Use Case 3

- [`salm_raw`](https://amrgen.org/reference/salm_raw.md) : Example
  Salmonella Genotype-Phenotype Data

- [`kleborate_raw`](https://amrgen.org/reference/kleborate_raw.md) :
  Example Kleborate Genotype Data from EuSCAPE project

- [`kleborate_raw_v313`](https://amrgen.org/reference/kleborate_raw_v313.md)
  : Example Kleborate v3.1.3 Genotype Data from EuSCAPE project

- [`kleborate_classes`](https://amrgen.org/reference/kleborate_classes.md)
  : Table mapping Kleborate drug class columns

- [`pheno_CLI_public`](https://amrgen.org/reference/pheno_CLI_public.md)
  :

  Clindamycin MIC data for 5914 *Staphylococcus aureus* isolates

- [`afp_CLI_public`](https://amrgen.org/reference/afp_CLI_public.md)
  : S. aureus Clindamycin Resistance Genotype Data

- [`ST_data_CLI`](https://amrgen.org/reference/ST_data_CLI.md) :

  ST data for *Staphylococcus aureus* genomes for clindamycin vignette

- [`DASSIM_pheno_raw`](https://amrgen.org/reference/DASSIM_pheno_raw.md)
  : DASSIM phenotype metadata

- [`btESBL_pheno`](https://amrgen.org/reference/btESBL_pheno.md) :
  blantyreESBL AST data

- [`DASSIM_geno`](https://amrgen.org/reference/DASSIM_geno.md) : DASSIM
  genotype data (AMRFinderPlus)

- [`NCBI_Ecoli_pheno_chl`](https://amrgen.org/reference/NCBI_Ecoli_pheno_chl.md)
  : NCBI AST for Escherichia coli tested against chloramphenicol

- [`MICROBIGGE_Ecoli_CHLR`](https://amrgen.org/reference/MICROBIGGE_Ecoli_CHLR.md)
  : MicroBIGG-E for E.coli containing all chloramphenicol resistance
  genes

- [`rgi_raw`](https://amrgen.org/reference/rgi_raw.md) : Example
  Resistance Gene Identifier (RGI) v6.0.6 Genotype Data

- [`rgi_short_name_table`](https://amrgen.org/reference/rgi_short_name_table.md)
  : Table mapping CARD/RGI Model ID and CARD Short Name

- [`rgi_drugs_table`](https://amrgen.org/reference/rgi_drugs_table.md) :
  Table mapping CARD/RGI drug class and antibiotic columns

- [`rgi_EuSCAPE_raw`](https://amrgen.org/reference/rgi_EuSCAPE_raw.md) :
  Example Resistance Gene Identifier (RGI) v6.0.6 Genotype Data from
  EuSCAPE project

- [`kp_mero_euscape`](https://amrgen.org/reference/kp_mero_euscape.md) :
  Meropenem Phenotype Data from EuSCAPE project

- [`kp_mero_amrfp`](https://amrgen.org/reference/kp_mero_amrfp.md) :
  Example AMRFinderPlus Genotype Data from EuSCAPE project

- [`sirscan_codes`](https://amrgen.org/reference/sirscan_codes.md) :
  SIRscan antibiotic code mapping
