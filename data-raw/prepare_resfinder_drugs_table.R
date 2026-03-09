# 1. Load your existing dataframe into an object EXACTLY named 'resfinder'
resfinder <- read.delim("resfinder_drug_classes_agents.tsv", check.names = FALSE)
# (Dataframe is in the structure: AFP_Subclass, ab, drug_class, drug_agent)

# 2. Use 'usethis' to save it as an .rda file in the data/ folder automatically
usethis::use_data(resfinder, overwrite = TRUE)
