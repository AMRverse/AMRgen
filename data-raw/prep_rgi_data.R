# Load example RGI results file
rgi_raw <- read_tsv("data-raw/rgi_v606_test.txt.gz")
usethis::use_data(rgi_raw, internal = FALSE, overwrite = TRUE)

# Load file mapping model ID column to CARD short names (https://card.mcmaster.ca/latest/data in aro_index.tsv)
rgi_short_name_table <- read_tsv("data-raw/rgi_modelid_shortname.tsv")
usethis::use_data(rgi_short_name_table, internal = FALSE, overwrite = TRUE)

# Load file mapping mapping CARD drug class / drug agents to standardised drug classes/names.
rgi_drugs_table <- read_tsv("data-raw/rgi_drug_classes_agents.tsv")
usethis::use_data(rgi_drugs_table, internal = FALSE, overwrite = TRUE)

# Load example RGI v6.0.6 results file using EuSCAPE (Klebsiella pneumoniae) data
rgi_EuSCAPE_raw <- read_tsv("data-raw/rgi_v606_EuSCAPE.txt.gz")
usethis::use_data(rgi_EuSCAPE_raw, internal = FALSE, overwrite = TRUE)

