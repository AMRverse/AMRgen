## General information
# Genotype files have been obtained using AMRfinderplus 4.2.7

####
# 2. DASSIM
####

## 2.5 Importing phenotype data

# Read in DASSIM metdata
DASSIM_pheno_raw <- read_tsv("data-raw/DASSIM_BSI_metadata.tsv")
usethis::use_data(DASSIM_pheno_raw, internal = FALSE, overwrite = TRUE)

# Read in the phenotype data from blantyre ESBL
btESBL_pheno <- readr::read_csv("data-raw/DASSIM_btESBL_AST.csv") %>% select(-row)
usethis::use_data(btESBL_pheno, internal = FALSE, overwrite = TRUE)


## 2.6 Importing genotype data

# Use import_amrfp() for AMRfinderplus data
DASSIM_geno <- import_amrfp("data-raw/DASSIM_BSI_Malawi_amrfinder_v4.0.23_geno.tsv")
usethis::use_data(DASSIM_geno, internal = FALSE, overwrite = TRUE)

####
# 3. NCBI
####

# 3.2  Reading in genotype and phenotype data

# import phenotype data
NCBI_Ecoli_pheno_chl <- import_pheno("data-raw/CHL_Ecoli_asts.tsv", format = "ncbi")
usethis::use_data(NCBI_Ecoli_pheno_chl, internal = FALSE, overwrite = TRUE)


# import genotype
MICROBIGGE_Ecoli_CHLR <- import_geno("data-raw/CHL-R_Ecoli_microbigge.tsv", format = "amrfp", sample_col = "BioSample")
usethis::use_data(MICROBIGGE_Ecoli_CHLR, internal = FALSE, overwrite = TRUE)
