# Get and Compare Antimicrobial Wild Type Distributions from EUCAST

These functions allow retrieval of antimicrobial wild type
distributions, live from [eucast.org](https://mic.eucast.org).

## Usage

``` r
get_eucast_amr_distribution(
  ab,
  mo = NULL,
  method = "MIC",
  as_freq_table = TRUE
)

get_eucast_mic_distribution(ab, mo = NULL, as_freq_table = TRUE)

get_eucast_disk_distribution(ab, mo = NULL, as_freq_table = TRUE)

compare_mic_with_eucast(mics, ab, mo = NULL)

compare_disk_with_eucast(disks, ab, mo = NULL)
```

## Arguments

- ab:

  antimicrobial, can be anything understood by
  [`ab_name()`](https://amr-for-r.org/reference/ab_property.html)

- mo:

  microorganism, can be anything understood by
  [`mo_name()`](https://amr-for-r.org/reference/mo_property.html) (can
  be left blank)

- method:

  either `"MIC"` or `"disk"`/`"diff"`

- as_freq_table:

  either `TRUE` (default) or `FALSE`, to return result as frequency
  table

- mics:

  MIC values, will be coerced with
  [`as.mic()`](https://amr-for-r.org/reference/as.mic.html)

- disks:

  Disk diffusion values, will be coerced with
  [`as.disk()`](https://amr-for-r.org/reference/as.disk.html)

## Details

The `compare_*_with_eucast()` functions allow to compare a user range
with EUCAST distributions. Use
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
on the output to visualise the results.

### Supported Antimicrobials

In December 2024, EUCAST had 176 distributions available, namely for
these antimicrobials:

Amikacin, amoxicillin, amoxicillin/clavulanic acid, amphotericin B,
ampicillin, ampicillin/sulbactam, anidulafungin, apramycin,
aspoxicillin, avilamycin, azithromycin, aztreonam, aztreonam/avibactam,
bacitracin, bedaquiline, benzylpenicillin, capreomycin, cefaclor,
cefadroxil, cefalexin, cefaloridine, cefalotin, cefapirin, cefazolin,
cefdinir, cefepime, cefepime/tazobactam, cefepime/zidebactam,
cefiderocol, cefixime, cefoperazone, cefoperazone/sulbactam, cefoselis,
cefotaxime, cefotetan, cefovecin, cefoxitin, cefpirome, cefpodoxime,
cefpodoxime/clavulanic acid, cefquinome, ceftaroline, ceftazidime,
ceftazidime/avibactam, ceftibuten, ceftiofur, ceftobiprole,
ceftolozane/tazobactam, ceftriaxone, cefuroxime, cephradine,
chloramphenicol, chlortetracycline, ciprofloxacin, clarithromycin,
clavulanic acid, clinafloxacin, clindamycin, clofazimine, cloxacillin,
colistin, cycloserine, dalbavancin, danofloxacin, daptomycin,
delafloxacin, delamanid, dicloxacillin, difloxacin, doripenem,
doxycycline, enrofloxacin, eravacycline, ertapenem, erythromycin,
ethambutol, ethionamide, faropenem, fidaxomicin, florfenicol,
flucloxacillin, fluconazole, flucytosine, flumequine, fosfomycin,
fusidic acid, gamithromycin, gatifloxacin, gemifloxacin, gentamicin,
imipenem, imipenem/relebactam, isavuconazole, isoniazid, itraconazole,
kanamycin, ketoconazole, lefamulin, levofloxacin, lincomycin, linezolid,
loracarbef, marbofloxacin, mecillinam, meropenem, meropenem/vaborbactam,
metronidazole, micafungin, minocycline, moxifloxacin, mupirocin,
nalidixic acid, narasin, neomycin, netilmicin, nitrofurantoin,
nitroxoline, norfloxacin, norvancomycin, ofloxacin, omadacycline,
orbifloxacin, oritavancin, oxacillin, oxolinic acid, oxytetracycline,
pefloxacin, phenoxymethylpenicillin, piperacillin,
piperacillin/tazobactam, pirlimycin, posaconazole, pradofloxacin,
pristinamycin, pyrazinamide, quinupristin/dalfopristin, retapamulin,
rezafungin, rifabutin, rifampicin, roxithromycin, secnidazole,
sitafloxacin, spectinomycin, spiramycin, streptomycin, sulbactam,
sulfadiazine, sulfamethoxazole, sulfisoxazole, tedizolid, teicoplanin,
telavancin, telithromycin, temocillin, terbinafine, tetracycline,
thiamphenicol, tiamulin, ticarcillin, ticarcillin/clavulanic acid,
tigecycline, tildipirosin, tilmicosin, tobramycin, trimethoprim,
trimethoprim/sulfamethoxazole, tulathromycin, tylosin, tylvalosin,
vancomycin, viomycin, and voriconazole.

For the current list, run
[`eucast_supported_ab_distributions()`](https://AMRverse.github.io/AMRgen/reference/eucast_supported_ab_distributions.md).

## Examples

``` r
get_eucast_mic_distribution("cipro")
#> # A tibble: 133 × 4
#>    microorganism                   microorganism_code   mic count
#>    <chr>                           <mo>               <mic> <int>
#>  1 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.002     0
#>  2 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.004     0
#>  3 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.008     0
#>  4 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.016     0
#>  5 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.030     0
#>  6 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.060     0
#>  7 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.125     0
#>  8 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.250    28
#>  9 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       0.500    27
#> 10 Actinobacillus pleuropneumoniae B_ACTNB_PLRP       1.000     1
#> # ℹ 123 more rows

# not returning as frequency table
get_eucast_mic_distribution("cipro", as_freq_table = FALSE)
#> # A tibble: 7 × 25
#>   microorganism microorganism_code `0.002` `0.004` `0.008` `0.016` `0.03` `0.06`
#>   <chr>         <mo>                 <int>   <int>   <int>   <int>  <int>  <int>
#> 1 Actinobacill… B_ACTNB_PLRP             0       0       0       0      0      0
#> 2 Mannheimia h… B_MNNHM_HMLY             0       0       0       0      0      0
#> 3 Pasteurella … B_PSTRL_MLTC             0       0       0       0      0      0
#> 4 Staphylococc… B_STPHY_AURS             0       0       0       0      0      0
#> 5 Staphylococc… B_STPHY_PSDN             0       0       0       0      1     16
#> 6 Streptococcu… B_STRPT_EQUI_EQUI        0       0       0       0      0      0
#> 7 Streptococcu… B_STRPT_EQUI_ZPDM        0       0       0       0      0      0
#> # ℹ 17 more variables: `0.125` <int>, `0.25` <int>, `0.5` <int>, `1` <int>,
#> #   `2` <int>, `4` <int>, `8` <int>, `16` <int>, `32` <int>, `64` <int>,
#> #   `128` <int>, `256` <int>, `512` <int>, distributions <int>,
#> #   observations <int>, ecoff <chr>, confidence_interval <chr>

# specify microorganism to only get results for that pathogen
get_eucast_mic_distribution("cipro", "K. pneumoniae")
#> # A tibble: 0 × 2
#> # ℹ 2 variables: mic <mic>, count <int>

get_eucast_disk_distribution("cipro", "K. pneumoniae")
#> # A tibble: 0 × 2
#> # ℹ 2 variables: disk_diffusion <dsk>, count <int>


# Plotting ----------------------------------------------------------------

mic_data <- get_eucast_mic_distribution("cipro", "K. pneumoniae")
mics <- rep(mic_data$mic, mic_data$count)
ggplot2::autoplot(mics, ab = "cipro", mo = "K. pneumoniae", title = "Look at my MICs!")
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Error in as.mic(): argument `x` must not be NULL


# Comparing With User Values ----------------------------------------------

my_mic_values <- AMR::random_mic(500)
comparison <- compare_mic_with_eucast(my_mic_values, ab = "cipro", mo = "K. pneumoniae")
#> Joining with `by = join_by(value)`
comparison
#> # A tibble: 20 × 3
#>    value   user eucast
#>  * <fct>  <int>  <int>
#>  1 0.0005    52      0
#>  2 0.001     40      0
#>  3 0.002     49      0
#>  4 0.004     29      0
#>  5 0.008     41      0
#>  6 0.016     25      0
#>  7 0.032     37      0
#>  8 0.064     28      0
#>  9 0.125     27      0
#> 10 0.25      27      0
#> 11 0.5       27      0
#> 12 1         16      0
#> 13 2         24      0
#> 14 4         20      0
#> 15 8         18      0
#> 16 16        11      0
#> 17 32        12      0
#> 18 64         7      0
#> 19 128        5      0
#> 20 >=256      5      0
#> Use ggplot2::autoplot() on this output to visualise.
ggplot2::autoplot(comparison)
#> Warning: Removed 20 rows containing missing values or values outside the scale range
#> (`geom_col()`).
```
