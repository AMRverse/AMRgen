# SIRscan antibiotic code mapping

Reference table mapping Bio-Rad SIRscan antibiotic codes to antibiotic
names recognised by
[`AMR::as.ab()`](https://amr-for-r.org/reference/as.ab.html). Used as
the default lookup in
[`import_sirscan_pheno()`](https://amrgen.org/reference/import_sirscan_pheno.md)
to translate codes that the AMR package does not recognise directly.
Users can supply their own table via the `sirscan_codes` argument if
their SIRscan configuration uses different code meanings.

## Usage

``` r
sirscan_codes
```

## Format

`sirscan_codes` A data frame with 11 rows and 2 columns:

Columns include:

- `sirscan_code`: The antibiotic code as it appears in SIRscan export
  files.

- `ab_name`: Antibiotic name string passed to
  [`AMR::as.ab()`](https://amr-for-r.org/reference/as.ab.html) for
  translation.

## Details

Confirmed codes (verified from Bio-Rad catalogue or collaborator
communication): `CF30` (cephalothin 30 mcg disk), `MA` (cefamandole),
`SSS` (sulfonamides), `SXT25` (trimethoprim/sulfamethoxazole 25 mcg
disk).

Codes mapped by French CASFM/SIRscan convention, not yet confirmed in
official Bio-Rad documentation: `S` (streptomycin), `K` (kanamycin), `C`
(chloramphenicol), `NALI` (nalidixic acid), `AZ` (azithromycin), `MM`
(minocycline), `ERM` (erythromycin).

Codes with no confirmed mapping (`CAN`, `ECM`, `ICM`) are absent from
this table; they will return `NA` from
[`AMR::as.ab()`](https://amr-for-r.org/reference/as.ab.html) with a
warning. Update this table or supply your own once the correct mappings
are known.
