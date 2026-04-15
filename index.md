# AMRgen

**AMRgen** is an open-source R package designed to **bridge the gap
between genotypic and phenotypic antimicrobial resistance (AMR) data**.
Developed as an extension to the [AMR R package](https://amr-for-r.org),
it provides tools to interpret AMR genes, integrate these findings with
antimicrobial susceptibility test (AST) data, and calculate
genotype-phenotype associations.

This package is developed in collaboration with the ESGEM-AMR Working
Group and is tailored for researchers and healthcare professionals
tackling AMR globally.

The [AMRgen website](https://amrgen.org) has full function
[documentation](https://amrgen.org/reference/index.html) and various
[vignettes](https://amrgen.org/articles/) working through analysing
geno/pheno data using key functions.

------------------------------------------------------------------------

## Key Features

- **Import Genotype and Phenotype Data**: Import from common formats
  (NCBI or EBI antibiogram format, VITEK, Sensititre, MicroScan,
  Phoenix, WHONet for phenotypes; AMRFinderPlus, ABRicate, Kleborate,
  RGI and AMRrules for genotypes.
- **Genotype-Phenotype Integration**: Links AMR gene presence with
  phenotypic resistance profiles, enabling deeper insights into
  resistance mechanisms.
- **Automated EUCAST MIC Distribution Integration**: Fetch MIC
  distribution data directly from [EUCAST](https://mic.eucast.org) for
  seamless comparison with local susceptibility data.
- **Visualisation**: Generate powerful UpSet plots to identify
  intersections of AMR gene presence and phenotypic resistance,
  highlighting multidrug resistance patterns.
- **Modular and Extensible**: Leverages the robust foundation of the AMR
  package, including antibiotic selectors and clinical breakpoint
  interpretations.
- **NCBI- and EBI-Compliant Export**: Export phenotype data to NCBI- and
  EBI-compliant antibiogram format.

------------------------------------------------------------------------

## Getting Started

To install and explore the package, follow the instructions below:

### Installation

It is best to restart R before running the installation to prevent
issues.

Install the latest version of this package with:

``` r
install.packages("remotes") # if you haven't already
remotes::install_github("AMRverse/AMRgen")
```

All required packages, including the [AMR
package](https://amr-for-r.org) if you don’t have it already, will be
installed automatically.

If you have issues, we recommend you install the latest version of the
AMR package directly, then try again to install AMRgen:

``` r
install.packages("remotes") # if you haven't already
remotes::install_github("msberends/AMR")
remotes::install_github("AMRverse/AMRgen")
```

If you still have trouble with installation please post an issue
[here](https://github.com/AMRverse/AMRgen/issues)! The package is new
and we want to make it as accessible as possible for new users.

## Quick Usage Examples

``` r
library(AMRgen)
```

### Investigate ciprofloxacin resistance vs quinolone genotype markers, via solo PPV and upset plots

Below is a quick example demonstrating the major functions with AMRgen
and how they can be used to compare ciprofloxacin resistance phenotypes
with quinolone genotype markers, in *E. coli*.

``` r
# Example public E. coli AST data from NCBI
#  (already imported via import_ncbi_pheno() and re-interpreted with as.sir())
ecoli_pheno

# Import matching E. coli AMRFinderPlus data from AllTheBacteria
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

# Calculate solo positive predictive value for ciprofloxacin resistance, for individual markers found solo
#  (for all quinolone-associated genotype markers)
soloPPV_cipro <- solo_ppv_analysis(ecoli_geno, ecoli_pheno, pheno_drug ="Ciprofloxacin", geno_class =c("Quinolones"), sir_col="pheno_clsi")

# Do upset plot of ciprofloxacin MIC vs quinolone genotype marker combinations
#  (for combinations observed at least 5 times)
cip_upset <- amr_upset(binary_matrix=soloPPV_cipro$amr_binary, min_set_size=5, assay="mic")

# Calculate positive predictive value for individual markers and combinations
cip_ppv <- ppv(binary_matrix=soloPPV_cipro$amr_binary, min_set_size=5)

# Do logistic regression of ciprofloxacin resistance as a function of presence/absence of quinolone-associated markers
#  (for markers observed at least 10 times)
models <- amr_logistic(binary_matrix=soloPPV_cipro$amr_binary, maf=10)

# Caclulate concordance of genotype and phenotype markers
eco_cip_matrix <- get_binary_matrix(ecoli_geno, ecoli_pheno, pheno_drug = "Ciprofloxacin", geno_class = "Quinolones", sir_col = "pheno_provided", keep_assay_values = TRUE, keep_assay_values_from = "mic")

concordance_cip <- concordance(eco_cip_matrix)
```

## Importing geno or pheno data

To learn how to download data from the public archives [NCBI or
EBI](http://amrgen.org/articles/DownloadGenoPhenoData.md)

You can see all the various phenotypic import functions for the
following formats by running:

``` r
?import_pheno
```

Currently, AMRgen supports the following phenotype formats:

- EBI
- NCBI
- VITEK
- MicroScan
- SensiTitre
- WHOnet

Phenotype data can also be exported in NCBI or EBI formats, for upload
to the public archives:

``` r
?export_ncbi_pheno
?export_ebi_pheno
```

## Example analyses using AMRgen

For a complete example of how to analyse your genotypic and phenotypic
data together, see [Analysing Geno-Pheno
Data](http://amrgen.org/articles/AnalysingGenoPhenoData.md).

Example analysis of a small *Salmonella enterica* dataset for
[ciprofloxacin
resistance](http://amrgen.org/articles/SalmonellaExamples.md)

Example analysis of clindamycin resistance in [*Staphylococcus
aureus*](http://amrgen.org/articles/StaphAureusClindamycin.md)

Example analysis using large-scale surviellance genotype and phenotype
data in [*Neiserria
gonohorroeae*](http://amrgen.org/articles/NeisseriaGonoExamples.md)

How to assess concordance of predicted phenotypes with [observed
phenotypes](http://amrgen.org/articles/Concordance.md)

For more see the various [vignettes](https://amrgen.org/articles/).

### Download reference MIC distribution from eucast.org and compare to example data

``` r
# Get MIC reference distribution for ciprofloxacin in E. coli
ecoli_cip_mic_data <- get_eucast_mic_distribution("cipro", "E. coli")

# Plot the reference distribution 
mics <- rep(ecoli_cip_mic_data$mic, ecoli_cip_mic_data$count)
ggplot2::autoplot(mics, ab = "cipro", mo = "E. coli", title = "E. coli cipro reference distribution")

# Compare reference distribution to example E. coli data
ecoli_cip <- ecoli_pheno$mic[ecoli_pheno$drug=="CIP"]
comparison <- compare_mic_with_eucast(ecoli_cip, ab = "cipro", mo = "E. coli")
comparison
ggplot2::autoplot(comparison)
```

## Contributions

Contributions are welcome! If you encounter issues or wish to suggest
new features, please open an issue or submit a pull request.

## Licence

This package is distributed under the GNU GPL-3.0 Licence. See `LICENSE`
for details.
