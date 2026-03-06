# N. gonorrhoeae Tetracycline Resistance Genotype Data Use Case 3

Raw concatenated output from AMRFinderPlus run on *Neisseria
gonorrhoeae* genome assemblies from 409 isolates collected in eastern
Spain (2021–2024). This object serves as the genotype input for
[import_amrfp](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
and is used to investigate the genetic basis of tetracycline resistance,
including the chromosomal *rpsJ* V57M mutation and the plasmid-borne
*tet(M)* gene.

## Usage

``` r
ngono_tet_geno_raw
```

## Format

`ngono_tet_geno_raw` A data frame with 4,428 rows and 24 columns:

- `Name`: Sample identifier (ENA run accession).

- `Protein id`: Protein identifier (all `NA`).

- `Contig id`: Assembly contig identifier.

- `Start`, `Stop`: Start and stop positions of the AMR element on the
  contig.

- `Strand`: Strand orientation (`+` or `-`).

- `Element symbol`: AMRFinderPlus element symbol (e.g. `rpsJ_V57M`,
  `tet(M)`).

- `Element name`: Full descriptive name of the AMR element.

- `Scope`: AMRFinderPlus scope (`core` or `plus`).

- `Type`: Type of AMR element (e.g. `AMR`).

- `Subtype`: Subtype of AMR element (e.g. `POINT`, `AMR`, `ALLELEX`).

- `Class`: Antibiotic class (e.g. `TETRACYCLINE`, `BETA-LACTAM`,
  `RIFAMYCIN`).

- `Subclass`: Antibiotic subclass (e.g. `TETRACYCLINE`, `CEPHALOSPORIN`,
  `RIFAMPIN`).

- `Method`: Detection method (e.g. `POINTX`, `BLASTX`, `ALLELEX`).

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

ENA project
[PRJEB83795](https://www.ebi.ac.uk/ena/browser/view/PRJEB83795). See
Sánchez-Serrano *et al.* (2026) <doi:10.1016/j.cmi.2025.12.026>.
