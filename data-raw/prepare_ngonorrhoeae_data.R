## General information
# Genotype files have been obtained using AMRfinderplus v4.0.23 with database version 2025-03-25.1

####
## Use case 1: Euro-GASP data
####

# Load raw genotype file
eurogasp_geno_raw <- read_tsv("data-raw/eurogasps_amrfp.tsv.gz")

# Load phenotype file
eurogasp_pheno_raw <- read_tsv("data-raw/eurogasps_MICs.tsv.gz")

usethis::use_data(eurogasp_geno_raw, internal = FALSE, overwrite = TRUE)
usethis::use_data(eurogasp_pheno_raw, internal = FALSE, overwrite = TRUE)

####
## Use case 2: PBP2 mutations
####

# Load raw genotype file
ngono_cro_geno_raw <- read_tsv("data-raw/ngono_cro_geno.tsv.gz")

# Load phenotype file
ngono_cro_pheno_raw <- read_tsv("data-raw/ngono_cro_pheno.tsv.gz")

# Convert to long format
usethis::use_data(ngono_cro_geno_raw, internal = FALSE, overwrite = TRUE)
usethis::use_data(ngono_cro_pheno_raw, internal = FALSE, overwrite = TRUE)

####
## Use case 3: Tetracycline resistance
####

# Load raw genotype file
ngono_tet_geno_raw <- read_tsv("data-raw/ngono_tet_geno.tsv.gz")

# Load phenotype file
ngono_tet_pheno_raw <- read_tsv("data-raw/ngono_tet_pheno.tsv.gz")

usethis::use_data(ngono_tet_geno_raw, internal = FALSE, overwrite = TRUE)
usethis::use_data(ngono_tet_pheno_raw, internal = FALSE, overwrite = TRUE)
