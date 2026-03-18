# S. aureus Clindamycin Resistance Genotype Data

Processed and filtered output from AMRFinderPlus run on *Staphylococcus
aureus* genome assemblies from AllTheBacteria. This object was processed
by [import_amrfp](https://amrgen.org/reference/import_amrfp.md) and is
used to investigate the genetic basis of clindamycin resistance.

## Usage

``` r
afp_CLI_public
```

## Format

`afp_CLI_public` A data frame with 43,287 rows and 42 columns:

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
  [import_amrfp](https://amrgen.org/reference/import_amrfp.md).

- `variation type`: Type of genetic variation interpreted by
  [import_amrfp](https://amrgen.org/reference/import_amrfp.md).

## Source

<https://github.com/AllTheBacteria/AllTheBacteria>
