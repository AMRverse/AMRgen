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

#' E. coli NCBI AST Example Data
#'
#' A subset of E. coli phenotype data from the NCBI AST browser.
#' @format A data frame with 10 rows and 21 columns representing unprocessed data from the NCBI AST browser.
#'
#' Columns include:
#' - `#BioSample`: Sample identifier.
#' - `Scientific name`: Species identifier.
#' - `Antibiotic`: Antibiotic name.
#' - `Testing standard`: Interpretation standard (EUCAST or CLSI).
#' - `Measurement sign`: Measurement sign (>, <, =, etc.) relating to MIC measurement.
#' - `MIC (mg/L)`: Minimum inhibitory concentration.
#' - `Disk diffusion (mm)`: Disk diffusion zone.
#' - `Resistance phenotype`: Resistance call (SIR) as submitted.
#' - ...: Additional metadata columns from the NCBI AST export.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"ecoli_ast_raw"


#' E. coli NCBI AST Example Data, Re-interpreted with AMR Package
#'
#' A subset of E. coli phenotype data from the NCBI AST browser.
#' @format A data frame with 4168 rows and 11 columns representing data from the NCBI AST browser, formatted and re-interpreted using [import_ast].
#'
#' Columns include:
#' - `id`: Sample identifier, imported from the `#BioSample` column in the raw input.
#' - `drug_agent`: Antibiotic code, interpreted from `Antibiotic` using `as.ab`, used to interpret `ecoff` and `pheno` columns.
#' - `mic`: Minimum inhibitory concentration, formatted using `as.mic`, used to interpret `ecoff` and `pheno` columns.
#' - `disk`: Disk diffusion zone, formatted using `as.disk`, used to interpret `ecoff` and `pheno` columns.
#' - `pheno_clsi`: S/I/R classification according to CLSI, interpreted using `as.sir`.
#' - `ecoff`: WT/NWT classification, interpreted using `as.sir`.
#' - `guideline`: Interpretation guidelines used to interpret `ecoff` and `pheno` columns.
#' - `method`: Test method, one of: `r toString(paste0('"', stats::na.omit(sort(unique(ecoli_ast$method))), "'"))`.
#' - `pheno_provided`: ??
#' - `spp_pheno`: Species identifier, interpreted from `Scientific name` using `as.mo`, used to interpret `ecoff` and `pheno` columns.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"ecoli_ast"


#' E. coli Genotype Example Data
#'
#' Genotypes called using AMRFinderPlus (v3.12.8, DB 2024-01-31.1), sourced from the AllTheBacteria project.
#' @format A data frame with 45228 rows and
#' 28 columns representing genotyping results
#' from AMRFinderPlus.
#'
#' Columns include:
#' - `Name`: Sample identifier.
#' - `Gene symbol`: Gene symbol in NCBI RefGene.
#' - `Hierarchy node`: Node in NCBI hierarchy.
#' - `Class`, `Subclass`: Drug class(es) associated with the marker (from NCBI RefGene).
#' - `% Coverage of reference sequence`, `% Identity to reference sequence`, `Accession of closest sequence`: Sequence match information.
#' - ...: Additional metadata columns from the AMRFinderPlus output.
#' @source <https://github.com/ncbi/amr/wiki/Interpreting-results>
"ecoli_geno_raw"


#' E. coli Ciprofloxacin MIC Distribution Example Data
#'
#' Ciprofloxacin MIC distributions for E. coli, calculated from public data
#' and compared with the EUCAST reference distribution.
#'
#' @format An object of class `compare_eucast` with 32 rows and
#' 3 columns. It provides MIC distributions
#' from EUCAST and public AST data extracted from [ecoli_ast]
#' in the form of counts per value.
#'
#' Columns include:
#' - `value`: MIC value.
#' - `user`: Count of samples with this MIC value, from the example data [ecoli_ast].
#' - `eucast`: Count of samples with this MIC value, downloaded from EUCAST (Feb 2026).
#' @source <https://mic.eucast.org/>
"ecoli_cip_vs_ref"


#' EUCAST Reference distribution for Ciprofloxacin in E. coli
#'
#' Data frame containing EUCAST reference distribution for ciprofloxacin in E. coli, downloaded using [get_eucast_mic_distribution].
#'
#' @format A data frame with 19 rows and 2 columns. It provides MIC distributions
#' from EUCAST in the form of counts per value.
#'
#' Columns include:
#' - `mic`: MIC value.
#' - `count`: Count of samples with this MIC value, downloaded from EUCAST (Feb 2026).
#' @source <https://mic.eucast.org/>
"ecoli_cip_mic_data"


#' S. aureus Example of Imported NCBI Phenotype Data
#'
#' Phenotypes sourced from NCBI Biosamples using the [download_ncbi_ast] function and imported to AMRgen phenotype table format.
#' @format `staph_ast_ncbi` A data frame with 143 rows and 19 columns representing all Staphylococcus aureus phenotyping results for amikacin and doxycycline downloaded from NCBI using [download_ncbi_ast], imported into AMRgen format using [import_ast].
#'
#' Columns include:
#' - `id`: Sample identifier.
#' - `drug_agent`: Antibiotic identifier, as class 'ab'.
#' - `mic`: MIC data, as class 'mic'.
#' - `disk`: Disk diffusion zone diameter data, as class 'disk'.
#' - `pheno_provided`, `pheno_eucast`: S/I/R phenotypes as downloaded from NCBI, and as re-interpreted from mic/disk measures against EUCAST 2024 breakpoints.
#' - ...: Additional data columns from NCBI.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"staph_ast_ncbi"


#' S. aureus Example of Raw Phenotype Data Downloaded from NCBI BioSamples via Entrez API
#'
#' Phenotypes sourced from NCBI Biosamples using the [download_ncbi_ast] function without reformating.
#' @format `staph_ast_ncbi_raw` A data frame with 143 rows and 13 columns representing all Staphylococcus aureus phenotyping results for amikacin and doxycycline.
#'
#' Columns include:
#' - `id`: Sample identifier.
#' - `Antibiotic`: Antibiotic name.
#' - `Resistance phenotype`: S/I/R phenotypes as downloaded from NCBI.
#' - ...: Additional data columns from NCBI.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"staph_ast_ncbi_raw"


#' S. aureus Example of Raw Phenotype Data Downloaded from NCBI via Google Cloud BigQuery
#'
#' Phenotypes sourced from NCBI via [query_ncbi_bq_ast] function, without reformating.
#' @format `staph_ast_ncbi_cloud_raw` A data frame with 142 rows and 11 columns representing all Staphylococcus aureus phenotyping results for amikacin and doxycycline.
#'
#' Columns include:
#' - `BioSample`: Sample identifier.
#' - `Antibiotic`: Antibiotic name.
#' - `Resistance phenotype`: S/I/R phenotypes as downloaded from NCBI.
#' - ...: Additional data columns from NCBI.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"staph_ast_ncbi_cloud_raw"


#' S. aureus Example of Raw Genotype Data Downloaded from NCBI via Google Cloud BigQuery
#'
#' AMRfinderplus genotypes sourced from NCBI via [query_ncbi_bq_ast] function, without reformating.
#' @format `staph_geno_ncbi_cloud_raw` A data frame with 4064 rows and 9 columns representing all Staphylococcus aureus genotyping results for markers associated with class aminoglycoside or tetracycline.
#'
#' Columns include:
#' - `biosample_acc`: Sample identifier.
#' - `scientific_name`: Organism name.
#' - `Gene symbol`, `Class`, `Subclass`, `Element type`, `Element subtype`, `Method`, `Hierarchy_node`: Key results fields from AMRfinderplus.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/microbigge/>
"staph_geno_ncbi_cloud_raw"


#' S. aureus Example of Imported EBI Phenotype Data
#'
#' Phenotypes sourced from EBI AMR Portal using the [download_ebi] function and imported to AMRgen phenotype table format.
#' @format `staph_ast_ebi` A data frame with 218 rows and 46 columns representing all Staphylococcus phenotyping results for amikacin and doxycycline downloaded from EBI using [download_ebi], and imported into AMRgen format using [import_ast].
#'
#' Columns include:
#' - `id`: Sample identifier.
#' - `drug_agent`: Antibiotic identifier, as class 'ab'.
#' - `mic`: MIC data, as class 'mic'.
#' - `disk`: Disk diffusion zone diameter data, as class 'disk'.
#' - `pheno_provided`, `pheno_eucast`, `pheno_clsi`, `ecoff`: S/I/R phenotypes as downloaded from EBI, and as re-interpreted from mic/disk measures against EUCAST 2024 breakpoints.
#' - ...: Additional data columns from EBI.
#' @source <https://www.ebi.ac.uk/amr>
"staph_ast_ebi"


#' S. aureus Example of Imported EBI Genotype Data
#'
#' Phenotypes sourced from EBI AMR Portal using the [download_ebi] function and imported to AMRgen phenotype table format.
#' @format `staph_geno_ebi` A data frame with 198344 rows and 34 columns representing all Staphylococcus genotyping results downloaded from EBI using [download_ebi], and imported into AMRgen format using [import_amrfp].
#'
#' Columns include:
#' - `id`: Sample identifier.
#' - `drug_agent`, `drug_class`: Antibiotic agent and class, determined by parsing AMRfinderplus `subclass` field in the downloaded file.
#' - `gene`, `node`, `marker`: gene symbol, parsed from `amr_element_symbol` field in the downloaded file.
#' - `mutation`: mutation within gene, parsed into HGVS nomenclature format from `amr_element_symbol` field in the downloaded file.
#' - `marker.label`: label for genotype marker, combining `gene` and `mutation` information (deletion variants represented as `"gene:-"`).
#' - ...: Additional data columns from EBI.
#' @source <https://www.ebi.ac.uk/amr>
"staph_geno_ebi"
