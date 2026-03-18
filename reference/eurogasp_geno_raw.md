# N. gonorrhoeae Euro-GASP Genotype Data Use Case 1

Raw concatenated output from AMRFinderPlus v4.0.23 (database
2025-03-25.1) run on *Neisseria gonorrhoeae* genome assemblies from
three Euro-GASP genomic surveys (2013, 2018, 2020). Assemblies were
generated with SPAdes v3.15.5 and quality-assessed with QUAST v5.1. This
object serves as the genotype input for
[import_amrfp](https://amrgen.org/reference/import_amrfp.md).

## Usage

``` r
eurogasp_geno_raw
```

## Format

`eurogasp_geno_raw` A data frame with 54,764 rows and 24 columns:

- `Name`: Sample identifier (ENA run accession).

- `Protein identifier`: Protein identifier (not used for *N.
  gonorrhoeae* assemblies; all `NA`).

- `Contig id`: Assembly contig identifier.

- `Start`, `Stop`: Start and stop positions of the AMR element on the
  contig.

- `Strand`: Strand orientation (`+` or `-`).

- `Element symbol`: AMRFinderPlus element symbol (gene name and/or
  mutation, e.g. `gyrA_S91F`).

- `Sequence name`: Full descriptive name of the AMR element.

- `Scope`: AMRFinderPlus scope (`core` or `plus`).

- `Element type`: Type of AMR element (e.g. `AMR`).

- `Element subtype`: Subtype of AMR element (e.g. `POINT`, `AMR`).

- `Class`: Antibiotic class (e.g. `QUINOLONE`, `BETA-LACTAM`).

- `Subclass`: Antibiotic subclass (e.g. `CEPHALOSPORIN`,
  `TETRACYCLINE`).

- `Method`: Detection method used by AMRFinderPlus (e.g. `POINTX`,
  `BLASTX`).

- `Target length`, `Reference sequence length`, `Alignment length`:
  Sequence length metrics (bp).

- `% Coverage of reference sequence`: Percentage of reference sequence
  covered by the alignment.

- `% Identity to reference sequence`: Percentage identity to the closest
  reference sequence.

- `Accession of closest sequence`: NCBI accession of the closest
  reference sequence.

- `Name of closest sequence`: Name of the closest reference sequence.

- `HMM id`, `HMM description`: HMM-based annotation fields (all `NA` for
  this dataset).

- `Hierarchy node`: Gene hierarchy node name, required for
  [import_amrfp](https://amrgen.org/reference/import_amrfp.md) (enabled
  by `--print_node` flag).

## Source

ENA projects
[PRJEB9227](https://www.ebi.ac.uk/ena/browser/view/PRJEB9227),
[PRJEB34068](https://www.ebi.ac.uk/ena/browser/view/PRJEB34068),
[PRJEB58139](https://www.ebi.ac.uk/ena/browser/view/PRJEB58139). See
also Harris *et al.* (2018)
<https://doi.org/10.1016/S1473-3099(18)30225-1>, Sánchez-Busó *et al.*
(2022) <https://doi.org/10.1016/S2666-5247(22)00044-1>, and Golparian
*et al.* (2024) <https://doi.org/10.1016/S2666-5247(23)00370-1>.
