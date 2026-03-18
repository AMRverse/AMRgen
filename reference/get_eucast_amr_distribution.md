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
[`eucast_supported_ab_distributions()`](https://amrgen.org/reference/eucast_supported_ab_distributions.md).

## Examples

``` r
if (FALSE) { # \dontrun{
get_eucast_mic_distribution("cipro")

# not returning as frequency table
get_eucast_mic_distribution("cipro", as_freq_table = FALSE)

# specify microorganism to only get results for that pathogen
get_eucast_mic_distribution("cipro", "K. pneumoniae")

get_eucast_disk_distribution("cipro", "K. pneumoniae")


# Plotting ----------------------------------------------------------------

mic_data <- get_eucast_mic_distribution("cipro", "K. pneumoniae")
mics <- rep(mic_data$mic, mic_data$count)
ggplot2::autoplot(mics, ab = "cipro", mo = "K. pneumoniae", title = "Look at my MICs!")


# Comparing With User Values ----------------------------------------------

my_mic_values <- AMR::random_mic(500)
comparison <- compare_mic_with_eucast(my_mic_values, ab = "cipro", mo = "K. pneumoniae")
comparison
ggplot2::autoplot(comparison)
} # }
```
