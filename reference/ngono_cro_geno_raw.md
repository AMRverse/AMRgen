# N. gonorrhoeae PBP2 Mutations Genotype Data Use Case 2

Raw concatenated output from AMRFinderPlus run on *Neisseria
gonorrhoeae* genome assemblies from multiple studies enriched for
ceftriaxone decreased susceptibility and resistance (total n = 2,101
genomes). This object serves as the genotype input for
[import_amrfp](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
and is used to investigate mosaic PBP2 (*penA*) variants associated with
ceftriaxone resistance.

## Usage

``` r
ngono_cro_geno_raw
```

## Format

`ngono_cro_geno_raw` A data frame with 23,885 rows and 24 columns:

- `Name`: Sample identifier (ENA run accession).

- `Protein id`: Protein identifier (all `NA`).

- `Contig id`: Assembly contig identifier.

- `Start`, `Stop`: Start and stop positions of the AMR element on the
  contig.

- `Strand`: Strand orientation (`+` or `-`).

- `Element symbol`: AMRFinderPlus element symbol (e.g. `pbp2`,
  `penA_A510V`).

- `Element name`: Full descriptive name of the AMR element.

- `Scope`: AMRFinderPlus scope (`core` or `plus`).

- `Element Type`: Type of AMR element (e.g. `AMR`).

- `Subtype`: Subtype of AMR element (e.g. `POINT`, `AMR`).

- `Class`: Antibiotic class (e.g. `BETA-LACTAM`).

- `Subclass`: Antibiotic subclass (e.g. `CEPHALOSPORIN`).

- `Method`: Detection method (e.g. `POINTX`, `BLASTX`).

- `Target length`, `Reference sequence length`, `Alignment length`:
  Sequence length metrics (bp).

- `% Coverage of reference`, `% Identity to reference`: Alignment
  quality metrics.

- `Closest reference accession`: NCBI accession of the closest reference
  sequence.

- `Closest reference name`: Name of the closest reference sequence.

- `HMM accession`, `HMM description`: HMM-based annotation fields (all
  `NA`).

- `Hierarchy node`: Gene hierarchy node name required for
  [import_amrfp](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md).

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
also Golparian *et al.* (2024)
<https://doi.org/10.1016/S2666-5247(23)00370-1>, Day *et al.* (2022)
<https://doi.org/10.2807/1560-7917.ES.2022.27.46.2200803>, Fifer *et
al.* (2024) <https://doi.org/10.1093/jac/dkae369>, van der Veen *et al.*
(2026) <https://doi.org/10.1093/cid/ciaf530>, Unemo *et al.* (2024)
<https://doi.org/10.1093/jac/dkae176>.
