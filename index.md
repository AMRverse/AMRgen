# AMRgen

**AMRgen** is an open-source R package designed to support systematic
AMR genotype-phenotype analysis. Developed as an extension to the [AMR R
package](https://amr-for-r.org), it provides functions to import and
harmonise genotypic data from common bioinformatics tools, alongside
phenotypic data from automated antimicrobial susceptibility testing
(AST) instruments and public repositories, and combine enotype-phenotype
data together in a single data structure. AMRgen supports common
analyses linking AST data to reference distributions, modelling
enotype-phenotype associations, quantifying concordance, and producing
publication-ready visualisations including UpSet plots that jointly
display genotypic marker combination frequencies and associated
phenotypic distributions.

This package is developed in collaboration with the [ESGEM-AMR Working
Group](https://esgem-amr.amrrules.org) and is tailored for researchers
and health professionals tackling AMR globally.

The [AMRgen website](https://amrgen.org) has full function
[documentation](https://amrgen.org/reference/index.html) and various
[vignettes](https://amrgen.org/articles/) working through analysing
geno/pheno data using key functions.

------------------------------------------------------------------------

## Key Features

- **Import genotype and phenotype data**: Import from common formats
  (NCBI or EBI antibiogram format, VITEK, Sensititre, MicroScan,
  Phoenix, WHONet for phenotypes; AMRFinderPlus, ABRicate, Kleborate,
  and CARD RGI for genotypes.
- **Fetch public genotype and phenotype data**: Download public data
  from NCBI or EBI.
- **Genotype-Phenotype Integration**: Links AMR gene presence with
  phenotypic resistance profiles, enabling deeper insights into
  resistance mechanisms.
- **Automated EUCAST reference distribution integration**: Fetch MIC and
  disk zone reference distribution data directly from
  [EUCAST](https://mic.eucast.org) for seamless comparison with local
  susceptibility data.
- **Visualisation**: Generate powerful UpSet plots to identify
  intersections of AMR gene presence and phenotypic resistance,
  highlighting multidrug resistance patterns.
- **Modular and extensible**: Leverages the robust foundation of the AMR
  package, including antibiotic selectors and clinical breakpoint
  interpretations.
- **NCBI- and EBI-Compliant Export**: Export phenotype data to
  antibiogram formats suitable for submission to NCBI or EBI along with
  genome data (linked by BioSample).

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
soloPPV_cipro <- solo_ppv(ecoli_geno, ecoli_pheno, pheno_drug ="Ciprofloxacin", geno_class =c("Quinolones"), sir_col="pheno_clsi")

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

Functions to import phenotypic data are documented in the
[import_pheno()](https://amrgen.org/reference/import_pheno.html)
function reference. Currently, AMRgen supports importing files in the
following formats:

- EBI
- NCBI
- VITEK
- MicroScan
- SensiTitre
- BD Phoenix
- WHOnet

The generic import function
[format_pheno()](https://amrgen.org/reference/format_pheno.html) is
provided to help import phenotype data from formats other than the
above. For a usage example see the vignettes on [Large-scale
regional/national surveillance
data](https://amrgen.org/articles/NeisseriaGonoExamples.html) and
[Custom stratification by isolate
source](https://amrgen.org/articles/SalmonellaExamples.html).

Functions to import genotypic data are documented in the
[import_geno()](https://amrgen.org/reference/import_geno.html) function
reference. Currently, AMRgen supports the following genotype data
formats:

- AMRFinderPlus
- Kleborate
- ABRicate (run with resfinder or ncbi databases)
- CARD RGI
- EBI AMR portal

## Fetching public geno or pheno data

To learn how to download data from the public archives (NCBI or EBI) see
the vignette: [Downloading Geno-Pheno
Data](http://amrgen.org/articles/DownloadGenoPhenoData.md)

Phenotype data can also be exported in NCBI or EBI formats, for upload
to the public archives, see documentation for
\[export_ebi_pheno()\](<https://amrgen.org/reference/export_ebi_pheno.html>\]
and
[export_ncbi_pheno()](https://amrgen.org/reference/export_ncbi_pheno.html)
functions.

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

## Vignettes: Example analyses using AMRgen

AMRgen includes a set of [vignettes](https://amrgen.org/articles/) to
illustrate package functionality.

For an overview of available functions, and how they can be used to
analyse genotypic and phenotypic data together, see \* [Analysing
Geno-Pheno Data](http://amrgen.org/articles/AnalysingGenoPhenoData.md).

Vignettes for specific tasks: \* [Downloading Geno-Pheno
Data](http://amrgen.org/articles/DownloadGenoPhenoData.md) from EBI or
NCBI databases \* [Assessing Geno-Pheno
Concordance](https://amrgen.org/articles/Concordance.html)

Vignettes using real-world examples to illustrate more complex
geno-pheno analyses: \* [Large-scale surveillance data for *Neiserria
gonohorroeae*](http://amrgen.org/articles/NeisseriaGonoExamples.md) \*
[Example with multiple *Salmonella enterica*
serovars](http://amrgen.org/articles/SalmonellaExamples.md) -
illustrating different ways to explore geno-pheno data by source,
genotypic marker count, etc \* [Analysing clindamycin resistance in
*Staphylococcus
aureus*](http://amrgen.org/articles/StaphAureusClindamycin.md) -
illustrating how to dig further into AMRFinderPlus genotype calls and
explore how different types of hits relate differently to phenotype \*
[Exploring catB3 deletion variants and impact on chloramphenicol
susceptibility in *Escherichia
coli*](https://amrgen.org/articles/DeletionVariantsCatB3.html) -
illustrating how to explore the impact of gene deletion variants on
phenotypes \* [Analysing meropenem resistance in *Klebsiella
pneumoniae*](https://amrgen.org/articles/ComparingGenotypers.html) -
exploring combined impacts of acquired genes and mutations, comparing
genotype calls from Kleborate, CARD RGI, AMRFinderPlus

## Contributions

Contributions are welcome! If you encounter issues or wish to suggest
new features, please open an issue or submit a pull request.

## Licence

This package is distributed under the GNU GPL-3.0 Licence. See `LICENSE`
for details.
