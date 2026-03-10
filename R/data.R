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


#' NCBI Subclass mapping to drug class
#'
#' Mapping of NCBI refgene / AMRfinderplus Subclass terms that are not present in the AMR package as drug class terms. Used internally when importing AMRfinderplus results into AMRgen genotype table format.
#' @format `amrfp_drugs_table` A data frame with 21 rows and 2 columns.
#'
#' Columns include:
#' - `AMRFP_Subclass`: NCBI term
#' - `drug_class`: Name of drug class as it should appear in imported genotype table.
"amrfp_drugs_table"

#' E. coli AST data from Mills et al 2022
#'
#' Phenotyping data published in Mills et al, Genome Medicine (2022) 14:147, downloaded via EBI AMR portal, and imported to AMRgen phenotype table format. Corresponding genotype data is available in `geno_eco_2075`.
#' @format `pheno_eco_2075` A data frame with 37350 rows and 37 columns representing MIC results for 2075 E. coli isolates tested against 18 drug agents, using BD Phoenix.
#'
#' Columns include:
#' - `id`: Sample identifier (BioSample).
#' - `drug_agent`: Antibiotic identifier, as class 'ab'.
#' - `mic`: MIC data, as class 'mic'.
#' - `pheno_provided`: S/I/R phenotypes as downloaded from EBI.
#' - ...: Additional data columns from EBI.
#' @source <https://www.ebi.ac.uk/amr>
"pheno_eco_2075"


#' E. coli genotype data from Mills et al 2022
#'
#' Genotyping data for isolates published in Mills et al, Genome Medicine (2022) 14:147, generated using AMRfinderplus v3.12.8 and downloaded from the [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria/tree/main/reproducibility/All-samples/AMR/AMRFinderPlus) project, and imported to AMRgen genotype table format. Corresponding MIC data is available in `pheno_eco_2075`.
#' @format `geno_eco_2075` A data frame with 56064 rows and 24 columns representing AMRfinderplus genotyping results for 2075 E. coli isolates.
#'
#' Columns include:
#' - `id`: Sample identifier.
#' - `drug_agent`, `drug_class`: Antibiotic agent and class, determined by parsing AMRfinderplus `subclass` field in the downloaded file.
#' - `gene`, `node`, `marker`: gene symbol, parsed from `amr_element_symbol` field in the downloaded file.
#' - `mutation`: mutation within gene, parsed into HGVS nomenclature format from `amr_element_symbol` field in the downloaded file.
#' - `marker.label`: label for genotype marker, combining `gene` and `mutation` information (deletion variants represented as `"gene:-"`).
#' - ...: Additional data columns from AMRfinderplus
#' @source <https://github.com/AllTheBacteria/AllTheBacteria>
"geno_eco_2075"


#' N. gonorrhoeae Euro-GASP Phenotype Data Use Case 1
#'
#' Minimum inhibitory concentration (MIC) data for *Neisseria gonorrhoeae* isolates
#' from three European Gonococcal Antimicrobial Surveillance Programme (Euro-GASP)
#' genomic surveys (2013, 2018, 2020), from
#' Harris *et al.* (2018) <doi:10.1016/S1473-3099(18)30225-1>,
#' Sánchez-Busó *et al.* (2022) <doi:10.1016/S2666-5247(22)00044-1>, and
#' Golparian *et al.* (2024) <doi:10.1016/S2666-5247(23)00370-1>.
#'
#' @format `eurogasp_pheno_raw` A data frame with 5,361 rows and 5 columns:
#' - `id`: Sample identifier (ENA run accession).
#' - `Azithromycin`: MIC value in mg/L for azithromycin (numeric; n = 5,055).
#' - `Ciprofloxacin`: MIC value in mg/L for ciprofloxacin (numeric; n = 5,360).
#' - `Cefixime`: MIC value in mg/L for cefixime (character, to preserve inequality prefixes such as `<0.016`).
#' - `Ceftriaxone`: MIC value in mg/L for ceftriaxone (numeric; n = 5,361).
#'
#' @source ENA projects [PRJEB9227](https://www.ebi.ac.uk/ena/browser/view/PRJEB9227),
#' [PRJEB34068](https://www.ebi.ac.uk/ena/browser/view/PRJEB34068),
#' [PRJEB58139](https://www.ebi.ac.uk/ena/browser/view/PRJEB58139).
"eurogasp_pheno_raw"


#' N. gonorrhoeae Euro-GASP Genotype Data Use Case 1
#'
#' Raw concatenated output from AMRFinderPlus v4.0.23 (database 2025-03-25.1) run on
#' *Neisseria gonorrhoeae* genome assemblies from three Euro-GASP genomic surveys
#' (2013, 2018, 2020). Assemblies were generated with SPAdes v3.15.5 and quality-assessed
#' with QUAST v5.1. This object serves as the genotype input for [import_amrfp].
#'
#' @format `eurogasp_geno_raw` A data frame with 54,764 rows and 24 columns:
#' - `Name`: Sample identifier (ENA run accession).
#' - `Protein identifier`: Protein identifier (not used for *N. gonorrhoeae* assemblies; all `NA`).
#' - `Contig id`: Assembly contig identifier.
#' - `Start`, `Stop`: Start and stop positions of the AMR element on the contig.
#' - `Strand`: Strand orientation (`+` or `-`).
#' - `Element symbol`: AMRFinderPlus element symbol (gene name and/or mutation, e.g. `gyrA_S91F`).
#' - `Sequence name`: Full descriptive name of the AMR element.
#' - `Scope`: AMRFinderPlus scope (`core` or `plus`).
#' - `Element type`: Type of AMR element (e.g. `AMR`).
#' - `Element subtype`: Subtype of AMR element (e.g. `POINT`, `AMR`).
#' - `Class`: Antibiotic class (e.g. `QUINOLONE`, `BETA-LACTAM`).
#' - `Subclass`: Antibiotic subclass (e.g. `CEPHALOSPORIN`, `TETRACYCLINE`).
#' - `Method`: Detection method used by AMRFinderPlus (e.g. `POINTX`, `BLASTX`).
#' - `Target length`, `Reference sequence length`, `Alignment length`: Sequence length metrics (bp).
#' - `% Coverage of reference sequence`: Percentage of reference sequence covered by the alignment.
#' - `% Identity to reference sequence`: Percentage identity to the closest reference sequence.
#' - `Accession of closest sequence`: NCBI accession of the closest reference sequence.
#' - `Name of closest sequence`: Name of the closest reference sequence.
#' - `HMM id`, `HMM description`: HMM-based annotation fields (all `NA` for this dataset).
#' - `Hierarchy node`: Gene hierarchy node name, required for [import_amrfp] (enabled by `--print_node` flag).
#'
#' @source ENA projects [PRJEB9227](https://www.ebi.ac.uk/ena/browser/view/PRJEB9227),
#' [PRJEB34068](https://www.ebi.ac.uk/ena/browser/view/PRJEB34068),
#' [PRJEB58139](https://www.ebi.ac.uk/ena/browser/view/PRJEB58139).
#' See also Harris *et al.* (2018) <doi:10.1016/S1473-3099(18)30225-1>,
#' Sánchez-Busó *et al.* (2022) <doi:10.1016/S2666-5247(22)00044-1>, and
#' Golparian *et al.* (2024) <doi:10.1016/S2666-5247(23)00370-1>.
"eurogasp_geno_raw"


#' N. gonorrhoeae PBP2 Mutations Phenotype Data Use Case 2
#'
#' Ceftriaxone MIC data for *Neisseria gonorrhoeae* isolates from multiple studies
#' enriched for decreased susceptibility and resistance to ceftriaxone, including
#' Euro-GASP 2020, ceftriaxone-resistant isolates from the United Kingdom and England,
#' isolates from Asia, and WHO reference genomes. Used to investigate mosaic PBP2
#' variants and their association with ceftriaxone MICs.
#'
#' @format `ngono_cro_pheno_raw` A data frame with 2,191 rows and 4 columns:
#' - `id`: Sample identifier (ENA run accession).
#' - `Ceftriaxone`: MIC value in mg/L for ceftriaxone (character, to preserve inequality prefixes).
#' - `penA`: *penA* mosaic allele type identifier (e.g. `60.001`, `237.001`).
#' - `study`: Source study identifier (e.g. `fifer2024`).
#'
#' @source ENA projects [PRJEB58139](https://www.ebi.ac.uk/ena/browser/view/PRJEB58139),
#' [PRJEB57389](https://www.ebi.ac.uk/ena/browser/view/PRJEB57389),
#' [PRJEB76977](https://www.ebi.ac.uk/ena/browser/view/PRJEB76977),
#' [PRJEB45627](https://www.ebi.ac.uk/ena/browser/view/PRJEB45627),
#' [PRJNA577446](https://www.ebi.ac.uk/ena/browser/view/PRJNA577446),
#' [PRJNA776899](https://www.ebi.ac.uk/ena/browser/view/PRJNA776899),
#' [PRJNA778600](https://www.ebi.ac.uk/ena/browser/view/PRJNA778600),
#' [PRJNA560592](https://www.ebi.ac.uk/ena/browser/view/PRJNA560592),
#' [PRJNA874857](https://www.ebi.ac.uk/ena/browser/view/PRJNA874857),
#' [PRJNA909328](https://www.ebi.ac.uk/ena/browser/view/PRJNA909328),
#' [PRJNA957547](https://www.ebi.ac.uk/ena/browser/view/PRJNA957547),
#' [PRJNA1161034](https://www.ebi.ac.uk/ena/browser/view/PRJNA1161034),
#' [PRJNA1189294](https://www.ebi.ac.uk/ena/browser/view/PRJNA1189294),
#' [PRJNA1067895](https://www.ebi.ac.uk/ena/browser/view/PRJNA1067895).
#' See also Golparian *et al.* (2024) <doi:10.1016/S2666-5247(23)00370-1>,
#' Day *et al.* (2022) <doi:10.2807/1560-7917.ES.2022.27.46.2200803>,
#' Fifer *et al.* (2024) <doi:10.1093/jac/dkae369>,
#' van der Veen *et al.* (2026) <doi:10.1093/cid/ciaf530>,
#' Unemo *et al.* (2024) <doi:10.1093/jac/dkae176>.
"ngono_cro_pheno_raw"

#' N. gonorrhoeae PBP2 Mutations Genotype Data Use Case 2
#'
#' Raw concatenated output from AMRFinderPlus run on *Neisseria gonorrhoeae* genome
#' assemblies from multiple studies enriched for ceftriaxone decreased susceptibility
#' and resistance (total n = 2,101 genomes). This object serves as the genotype input
#' for [import_amrfp] and is used to investigate mosaic PBP2 (*penA*) variants
#' associated with ceftriaxone resistance.
#'
#' @format `ngono_cro_geno_raw` A data frame with 23,885 rows and 24 columns:
#' - `Name`: Sample identifier (ENA run accession).
#' - `Protein id`: Protein identifier (all `NA`).
#' - `Contig id`: Assembly contig identifier.
#' - `Start`, `Stop`: Start and stop positions of the AMR element on the contig.
#' - `Strand`: Strand orientation (`+` or `-`).
#' - `Element symbol`: AMRFinderPlus element symbol (e.g. `pbp2`, `penA_A510V`).
#' - `Element name`: Full descriptive name of the AMR element.
#' - `Scope`: AMRFinderPlus scope (`core` or `plus`).
#' - `Element Type`: Type of AMR element (e.g. `AMR`).
#' - `Subtype`: Subtype of AMR element (e.g. `POINT`, `AMR`).
#' - `Class`: Antibiotic class (e.g. `BETA-LACTAM`).
#' - `Subclass`: Antibiotic subclass (e.g. `CEPHALOSPORIN`).
#' - `Method`: Detection method (e.g. `POINTX`, `BLASTX`).
#' - `Target length`, `Reference sequence length`, `Alignment length`: Sequence length metrics (bp).
#' - `% Coverage of reference`, `% Identity to reference`: Alignment quality metrics.
#' - `Closest reference accession`: NCBI accession of the closest reference sequence.
#' - `Closest reference name`: Name of the closest reference sequence.
#' - `HMM accession`, `HMM description`: HMM-based annotation fields (all `NA`).
#' - `Hierarchy node`: Gene hierarchy node name required for [import_amrfp].
#'
#' @source ENA projects
#' [PRJEB58139](https://www.ebi.ac.uk/ena/browser/view/PRJEB58139),
#' [PRJEB57389](https://www.ebi.ac.uk/ena/browser/view/PRJEB57389),
#' [PRJEB76977](https://www.ebi.ac.uk/ena/browser/view/PRJEB76977),
#' [PRJEB45627](https://www.ebi.ac.uk/ena/browser/view/PRJEB45627),
#' [PRJNA577446](https://www.ebi.ac.uk/ena/browser/view/PRJNA577446),
#' [PRJNA776899](https://www.ebi.ac.uk/ena/browser/view/PRJNA776899),
#' [PRJNA778600](https://www.ebi.ac.uk/ena/browser/view/PRJNA778600),
#' [PRJNA560592](https://www.ebi.ac.uk/ena/browser/view/PRJNA560592),
#' [PRJNA874857](https://www.ebi.ac.uk/ena/browser/view/PRJNA874857),
#' [PRJNA909328](https://www.ebi.ac.uk/ena/browser/view/PRJNA909328),
#' [PRJNA957547](https://www.ebi.ac.uk/ena/browser/view/PRJNA957547),
#' [PRJNA1161034](https://www.ebi.ac.uk/ena/browser/view/PRJNA1161034),
#' [PRJNA1189294](https://www.ebi.ac.uk/ena/browser/view/PRJNA1189294),
#' [PRJNA1067895](https://www.ebi.ac.uk/ena/browser/view/PRJNA1067895).
#' See also Golparian *et al.* (2024) <doi:10.1016/S2666-5247(23)00370-1>,
#' Day *et al.* (2022) <doi:10.2807/1560-7917.ES.2022.27.46.2200803>,
#' Fifer *et al.* (2024) <doi:10.1093/jac/dkae369>,
#' van der Veen *et al.* (2026) <doi:10.1093/cid/ciaf530>,
#' Unemo *et al.* (2024) <doi:10.1093/jac/dkae176>.
"ngono_cro_geno_raw"


#' N. gonorrhoeae Tetracycline Resistance Phenotype Data Use Case 3
#'
#' Tetracycline MIC data for 409 *Neisseria gonorrhoeae* isolates collected in
#' Eastern Spain between 2021 and 2024, used to investigate genetic determinants
#' of tetracycline resistance and implications for doxy-PEP strategies.
#' Described in Sánchez-Serrano *et al.* (2026) <doi:10.1016/j.cmi.2025.12.026>.
#'
#' @format `ngono_tet_pheno_raw` A data frame with 409 rows and 2 columns:
#' - `id`: Sample identifier (ENA run accession).
#' - `Tetracycline`: MIC value in mg/L for tetracycline (numeric).
#'
#' @source ENA project [PRJEB83795](https://www.ebi.ac.uk/ena/browser/view/PRJEB83795).
#' See Sánchez-Serrano *et al.* (2026) <doi:10.1016/j.cmi.2025.12.026>.
"ngono_tet_pheno_raw"


#' N. gonorrhoeae Tetracycline Resistance Genotype Data Use Case 3
#'
#' Raw concatenated output from AMRFinderPlus run on *Neisseria gonorrhoeae* genome
#' assemblies from 409 isolates collected in eastern Spain (2021–2024). This object
#' serves as the genotype input for [import_amrfp] and is used to investigate
#' the genetic basis of tetracycline resistance, including the chromosomal *rpsJ* V57M
#' mutation and the plasmid-borne *tet(M)* gene.
#'
#' @format `ngono_tet_geno_raw` A data frame with 4,428 rows and 24 columns:
#' - `Name`: Sample identifier (ENA run accession).
#' - `Protein id`: Protein identifier (all `NA`).
#' - `Contig id`: Assembly contig identifier.
#' - `Start`, `Stop`: Start and stop positions of the AMR element on the contig.
#' - `Strand`: Strand orientation (`+` or `-`).
#' - `Element symbol`: AMRFinderPlus element symbol (e.g. `rpsJ_V57M`, `tet(M)`).
#' - `Element name`: Full descriptive name of the AMR element.
#' - `Scope`: AMRFinderPlus scope (`core` or `plus`).
#' - `Type`: Type of AMR element (e.g. `AMR`).
#' - `Subtype`: Subtype of AMR element (e.g. `POINT`, `AMR`, `ALLELEX`).
#' - `Class`: Antibiotic class (e.g. `TETRACYCLINE`, `BETA-LACTAM`, `RIFAMYCIN`).
#' - `Subclass`: Antibiotic subclass (e.g. `TETRACYCLINE`, `CEPHALOSPORIN`, `RIFAMPIN`).
#' - `Method`: Detection method (e.g. `POINTX`, `BLASTX`, `ALLELEX`).
#' - `Target length`, `Reference sequence length`, `Alignment length`: Sequence length metrics (bp).
#' - `% Coverage of reference`, `% Identity to reference`: Alignment quality metrics.
#' - `Closest reference accession`: NCBI accession of the closest reference sequence.
#' - `Closest reference name`: Name of the closest reference sequence.
#' - `HMM accession`, `HMM description`: HMM-based annotation fields (all `NA`).
#' - `Hierarchy node`: Gene hierarchy node name required for [import_amrfp].
#'
#' @source ENA project [PRJEB83795](https://www.ebi.ac.uk/ena/browser/view/PRJEB83795).
#' See Sánchez-Serrano *et al.* (2026) <doi:10.1016/j.cmi.2025.12.026>.
"ngono_tet_geno_raw"

#' Example Salmonella Genotype-Phenotype Data
#'
#' Raw genotype-phenotype data for *Salmonella enterica* genomes, one row per sample.
#'
#' @format `salm_raw` A data frame with 115 rows and 7 columns:
#' - `Sample`: Sample identifier
#' - `Source`, `Serovar`: Non-AMR related information about each isolate
#' - `CpL_Genotype`: List of genotypic parkers (separated by `;`)
#' - `Ciprofloxacin`, `Levofloxacin`, `Moxifloxacin`: MIC value for each drug
"salm_raw"

#' Example Kleborate Genotype Data
#'
#' Raw Kleborate results file for *Klebsiella pneumoniae* genomes, one row per sample.
#'
#' @format `kleborate_raw` A data frame with 1,689 rows and 122 columns:
#' - `strain`: Sample identifier
#' - ...: Kleborate results columns
"kleborate_raw"
