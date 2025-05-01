# AMRgen

R Package for Genetic and Phenotypic Resistance Interpretation 

## Overview

**AMRgen** is an open-source R package designed to bridge the gap between genotypic and phenotypic antimicrobial resistance (AMR) data. Developed as an extension to the [AMR R package](https://github.com/msberends/AMR), it provides tools to interpret AMR genes, integrate these findings with antimicrobial susceptibility test (AST) data, and calculate genotype-phenotype associations.

This package is developed in collaboration with the ESGEM-AMR Working Group and is tailored for researchers and healthcare professionals tackling AMR globally.

The [AMRgen website](https://interpretamr.github.io/AMRgen/index.html) has full function [documentation](https://interpretamr.github.io/AMRgen/reference/index.html) and a [Vignette](https://interpretamr.github.io/AMRgen/articles/AnalysingGenoPhenoData.html) working through analysing geno/pheno data using key functions.

------------------------------------------------------------------------

## Key Features

-   **Import Genotype and Phenotype Data**: Import from common formats (NCBI AST for phenotypes; AMRfinderplus and hAMRonization for genotypes
-   **Genotype-Phenotype Integration**: Links AMR gene presence with phenotypic resistance profiles, enabling deeper insights into resistance mechanisms.
-   **Automated EUCAST MIC Distribution Integration**: Fetch MIC distribution data directly from [EUCAST](https://mic.eucast.org) for seamless comparison with local susceptibility data.
-   **Visualisation**: Generate powerful UpSet plots to identify intersections of AMR gene presence and phenotypic resistance, highlighting multidrug resistance patterns.
-   **Modular and Extensible**: Leverages the robust foundation of the AMR package, including antibiotic selectors and clinical breakpoint interpretations.

> Planned for development
-   **NCBI-Compliant Export**: Export phenotype data to NCBI-compliant antibiogram format.
-   **Expanded Data Import**: Import and parse phenotype data from other tools (e.g. CARD, ResFinder).

------------------------------------------------------------------------

## Getting Started

To install and explore the package, follow the instructions below:

### Installation
Note that this package requires the latest version of the `AMR` package (still in beta).

Install the latest version of the `AMR` package with:
```r
install.packages("remotes") # if you haven't already
remotes::install_github("msberends/AMR")
```

Then install this package
```r
# Install from GitHub
remotes::install_github("interpretAMR/AMRgen")
```

It is best to restart R before running the installation. If you didn't do this and/or you encounter issues with the examples below after install, it may help to also restart after the install and start fresh with the examples below.

## Quick Usage Examples

```r
library(AMRgen)
```

### Investigate ciprofloxacin resistance vs quinolone genotype markers, via solo PPV and upset plots

```r
# Example public E. coli AST data from NCBI
#  (already imported via import_ncbi_ast() and re-interpreted with as.sir())
ecoli_ast

# Import matching E. coli AMRfinderplus data from AllTheBacteria
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

# Calculate positive predictive value for ciprofloxacin resistance
#  (for all quinolone-associated genotype markers)
soloPPV_cipro <- solo_ppv_analysis(ecoli_geno, ecoli_ast, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="pheno")

# Do upset plot of ciprofloxacin MIC vs quinolone genotype marker combinations
#  (for combinations observed at least 5 times)
cip_upset <- amr_upset(soloPPV_cipro$amr_binary, min_set_size=5, assay="mic", order="value")

# Do logistic regression of ciprofloxacin resistance as a function of presence/absence of quinolone-associated markers
#  (for markers observed at least 10 times)
models <- amr_logistic(geno_table = import_amrfp(ecoli_geno_raw, "Name"), pheno_table = ecoli_ast, 
                       antibiotic = "Ciprofloxacin", drug_class_list = c("Quinolones"), maf=10)

```

### Download reference MIC distribution from eucast.org and compare to example data

```r
# Get MIC reference distribution for ciprofloxacin in E. coli
ecoli_cip_mic_data <- get_eucast_mic_distribution("cipro", "E. coli")

# Plot the reference distribution 
mics <- rep(ecoli_cip_mic_data$mic, ecoli_cip_mic_data$count)
ggplot2::autoplot(mics, ab = "cipro", mo = "E. coli", title = "E. coli cipro reference distribution")

# Compare reference distribution to example E. coli data
ecoli_cip <- ecoli_ast$mic[ecoli_ast$drug_agent=="CIP"]
comparison <- compare_mic_with_eucast(ecoli_cip, ab = "cipro", mo = "E. coli")
comparison
ggplot2::autoplot(comparison)
```

For more see the [Vignette](https://interpretamr.github.io/AMRgen/articles/AnalysingGenoPhenoData.html).


## Contributions

Contributions are welcome! If you encounter issues or wish to suggest new features, please open an issue or submit a pull request.

## Licence

This package is distributed under the GNU GPL-3.0 Licence. See `LICENSE` for details.
