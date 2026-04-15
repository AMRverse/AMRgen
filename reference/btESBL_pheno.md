# blantyreESBL AST data

The antimicrobial susceptibility testing (AST) data from the
blantyreESBL Github from Dr. Joseph Lewis. Antimicrobial sensitivity
testing (AST) was carried out on a subset of isolates using the
disc-diffusion method using British Society for Antimicrobial
Chemotherapy (BSAC) guidelines (https://bsac.org.uk/). AST was carried
out for meropenem, amikacin, chloramphenicol, ciprofloxacin,
co-trimoxazole and gentamicin. However this dataset only contains the
raw phenotype data (S/R).

## Usage

``` r
btESBL_pheno
```

## Format

`btESBL_pheno` A data frame with 609 rows and 9 columns:

- `...1`: Row count.

- `supplier_name`: Strain ID.

- `organism`: organism.

- `amikacin`: amikacin antimicrobial susceptibility phenotype (S/R)

- `chloramphenicol`: chloramphenicol antimicrobial susceptibility
  phenotype (S/R)

- `ciprofloxacin`: ciprofloxacin antimicrobial susceptibility phenotype
  (S/R)

- `cotrimoxazole`: cotrimoxazole antimicrobial susceptibility phenotype
  (S/R)

- `gentamicin`: gentamicin antimicrobial susceptibility phenotype (S/R)

- `meropenem`: meropenem antimicrobial susceptibility phenotype (S/R)

## Source

<https://github.com/joelewis101/blantyreESBL/raw/refs/heads/main/data/btESBL_pheno.rda>
