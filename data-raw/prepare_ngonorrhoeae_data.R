## General information
# Genotype files have been obtained using AMRFinderPlus v4.0.23 with database version 2025-03-25.1

####
## Use case 1: Euro-GASP data
####

# Load raw genotype file
eurogasp_geno_raw <- read_tsv("data-raw/eurogasps_amrfp.tsv.gz")

# Load phenotype file
eurogasp_pheno_raw <- read_tsv("data-raw/eurogasps_MICs.tsv.gz")

usethis::use_data(eurogasp_geno_raw, internal = FALSE, overwrite = TRUE)
usethis::use_data(eurogasp_pheno_raw, internal = FALSE, overwrite = TRUE)

# Load comparison of azithromycin data to EUCAST reference distribution
# (retrieved 15/3/2026)
azm_comparison <- read_tsv("data-raw/Ngono_azm_comparison.tsv")
azm_comparison <- structure(azm_comparison, class = c("compare_eucast", class(azm_comparison)))
usethis::use_data(azm_comparison, internal = FALSE, overwrite = TRUE)

# Load comparison of ciprofloxacin data to EUCAST reference distribution
# (retrieved 15/3/2026)
cip_comparison <- read_tsv("data-raw/Ngono_cip_comparison.tsv")
cip_comparison <- structure(cip_comparison, class = c("compare_eucast", class(cip_comparison)))
usethis::use_data(cip_comparison, internal = FALSE, overwrite = TRUE)

# Load comparison of ceftriaxone data to EUCAST reference distribution
# (retrieved 15/3/2026)
cro_comparison <- read_tsv("data-raw/Ngono_cro_comparison.tsv")
cro_comparison <- structure(cro_comparison, class = c("compare_eucast", class(cro_comparison)))
usethis::use_data(cro_comparison, internal = FALSE, overwrite = TRUE)

# Load comparison of cefixime data to EUCAST reference distribution
# (retrieved 15/3/2026)
cfm_comparison <- read_tsv("data-raw/Ngono_cfm_comparison.tsv")
cfm_comparison <- structure(cfm_comparison, class = c("compare_eucast", class(cfm_comparison)))
usethis::use_data(cfm_comparison, internal = FALSE, overwrite = TRUE)

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
