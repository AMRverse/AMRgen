# Load example kleborate results file
kleborate_raw <- read_tsv("data-raw/EuSCAPE_kleborate.tsv.gz")
usethis::use_data(kleborate_raw, internal = FALSE, overwrite = TRUE)

# load file mapping kleborate AMR columns to AMR pkg drug classes
kleborate_classes <- read_tsv("data-raw/kleborate_drug_classes.tsv")
usethis::use_data(kleborate_classes, internal = FALSE, overwrite = TRUE)
