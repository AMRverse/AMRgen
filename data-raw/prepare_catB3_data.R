## General information
# Genotype files have been obtained using AMRfinderplus 4.2.7 

####
# 2. DASSIM
####

## 2.5 Importing phenotype data

# Read in DASSIM metdata 
DASSIM_pheno = read_tsv("data-raw/DASSIM_BSI_metadata.tsv")
usethis::use_data(DASSIM_pheno, internal = FALSE, overwrite = TRUE)

# Read in the AST data from blantyre ESBL
btESBL_AST = readr::read_csv("data-raw/DASSIM_btESBL_AST.csv")
usethis::use_data(btESBL_AST, internal = FALSE, overwrite = TRUE)


## 2.6 Importing genotype data

# Use import_amrfp() for AMRfinderplus data
DASSIM_geno = import_amrfp("data-raw/DASSIM_BSI_Malawi_amrfinder_v4.0.23_geno.tsv")
usethis::use_data(DASSIM_geno, internal = FALSE, overwrite = TRUE)

####
# 3. NCBI
####

# 3.2  Reading in genotype and phenotype data

# import phenotype data
NCBI_AST_CHL = read_tsv("data-raw/CHL_Ecoli_asts.tsv")
NCBI_AST_CHL = NCBI_AST_CHL %>% rename("BioSample" = "#BioSample")
usethis::use_data(NCBI_AST_CHL, internal = FALSE, overwrite = TRUE)


# import genotype
MICROBIGGE_CATB3 = read_tsv("data-raw/catB3_Ecoli_microbigge.tsv")
usethis::use_data(MICROBIGGE_CATB3, internal = FALSE, overwrite = TRUE)
MICROBIGGE_CHLR = read_tsv("data-raw/CHL-R_Ecoli_microbigge.tsv")
usethis::use_data(MICROBIGGE_CHLR, internal = FALSE, overwrite = TRUE)
