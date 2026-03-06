# N. gonorrhoeae PBP2 Mutations Phenotype Data Use Case 2

Ceftriaxone MIC data for *Neisseria gonorrhoeae* isolates from multiple
studies enriched for decreased susceptibility and resistance to
ceftriaxone, including Euro-GASP 2020, ceftriaxone-resistant isolates
from the United Kingdom and England, isolates from Asia, and WHO
reference genomes. Used to investigate mosaic PBP2 variants and their
association with ceftriaxone MICs.

## Usage

``` r
ngono_cro_pheno_raw
```

## Format

`ngono_cro_pheno_raw` A data frame with 2,191 rows and 4 columns:

- `id`: Sample identifier (ENA run accession).

- `Ceftriaxone`: MIC value in mg/L for ceftriaxone (character, to
  preserve inequality prefixes).

- `penA`: *penA* mosaic allele type identifier (e.g. `60.001`,
  `237.001`).

- `study`: Source study identifier (e.g. `fifer2024`).

## Source

ENA projects
[PRJEB58139](https://www.ebi.ac.uk/ena/browser/view/PRJEB58139),
[PRJEB57389](https://www.ebi.ac.uk/ena/browser/view/PRJEB57389),
[PRJEB76977](https://www.ebi.ac.uk/ena/browser/view/PRJEB76977),
[PRJEB45627](https://www.ebi.ac.uk/ena/browser/view/PRJEB45627),
[PRJNA577446](https://www.ebi.ac.uk/ena/browser/view/PRJNA577446),
[PRJNA776899](https://www.ebi.ac.uk/ena/browser/view/PRJNA776899),
[PRJNA778600](https://www.ebi.ac.uk/ena/browser/view/PRJNA778600),
[PRJNA560592](https://www.ebi.ac.uk/ena/browser/view/PRJNA560592),
[PRJNA874857](https://www.ebi.ac.uk/ena/browser/view/PRJNA874857),
[PRJNA909328](https://www.ebi.ac.uk/ena/browser/view/PRJNA909328),
[PRJNA957547](https://www.ebi.ac.uk/ena/browser/view/PRJNA957547),
[PRJNA1161034](https://www.ebi.ac.uk/ena/browser/view/PRJNA1161034),
[PRJNA1189294](https://www.ebi.ac.uk/ena/browser/view/PRJNA1189294),
[PRJNA1067895](https://www.ebi.ac.uk/ena/browser/view/PRJNA1067895). See
also Golparian *et al.* (2024) <doi:10.1016/S2666-5247(23)00370-1>, Day
*et al.* (2022) <doi:10.2807/1560-7917.ES.2022.27.46.2200803>, Fifer *et
al.* (2024) <doi:10.1093/jac/dkae369>, van der Veen *et al.* (2026)
<doi:10.1093/cid/ciaf530>, Unemo *et al.* (2024)
<doi:10.1093/jac/dkae176>.
