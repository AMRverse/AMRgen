# Retrieve Available Antimicrobial Wild Type Distributions from EUCAST

Run this function to get an updated list of antimicrobial distributions
currently supported by EUCAST. This retrieves live info from
<https://mic.eucast.org>.

## Usage

``` r
eucast_supported_ab_distributions(...)
```

## Arguments

- ...:

  Arguments passed on to the function, currently unused.

## Examples

``` r
eucast_supported_ab_distributions()
#>                             AMK                             AMX 
#>                      "Amikacin"                   "Amoxicillin" 
#>                             AMC                             AMB 
#>   "Amoxicillin/clavulanic acid"                "Amphotericin B" 
#>                             AMP                             SAM 
#>                    "Ampicillin"          "Ampicillin/sulbactam" 
#>                             SAM                             APR 
#>          "Ampicillin/sulbactam"                     "Apramycin" 
#>                             APX                             AVI 
#>                  "Aspoxicillin"                    "Avilamycin" 
#>                             AZM                             ATM 
#>                  "Azithromycin"                     "Aztreonam" 
#>                             AZA                             BAC 
#>           "Aztreonam/avibactam"                    "Bacitracin" 
#>                             BDQ                             PEN 
#>                   "Bedaquiline"              "Benzylpenicillin" 
#>                             CAP                             CEC 
#>                   "Capreomycin"                      "Cefaclor" 
#>                             CFR                             LEX 
#>                    "Cefadroxil"                     "Cefalexin" 
#>                             RID                             CEP 
#>                  "Cefaloridine"                     "Cefalotin" 
#>                             HAP                             CZO 
#>                     "Cefapirin"                     "Cefazolin" 
#>                             CDR                             FEP 
#>                      "Cefdinir"                      "Cefepime" 
#>                             FPE                             FPT 
#>       "Cefepime/enmetazobactam"           "Cefepime/tazobactam" 
#>                             FPZ                             FDC 
#>           "Cefepime/zidebactam"                   "Cefiderocol" 
#>                             CFM                             CFP 
#>                      "Cefixime"                  "Cefoperazone" 
#>                             CSL                             CSE 
#>        "Cefoperazone/sulbactam"                     "Cefoselis" 
#>                             CTX                             CTT 
#>                    "Cefotaxime"                     "Cefotetan" 
#>                             FOV                             FOX 
#>                     "Cefovecin"                     "Cefoxitin" 
#>                             CPO                             CPD 
#>                     "Cefpirome"                   "Cefpodoxime" 
#>                             CDC                             CEQ 
#>   "Cefpodoxime/clavulanic acid"                    "Cefquinome" 
#>                             CPT                             CAZ 
#>                   "Ceftaroline"                   "Ceftazidime" 
#>                             CZA                             CTB 
#>         "Ceftazidime/avibactam"                    "Ceftibuten" 
#>                             TIO                             BPR 
#>                     "Ceftiofur"                  "Ceftobiprole" 
#>                             CZT                             CZT 
#>        "Ceftolozane/tazobactam"        "Ceftolozane/tazobactam" 
#>                             CRO                             CXM 
#>                   "Ceftriaxone"                    "Cefuroxime" 
#>                             CED                             CHL 
#>                    "Cephradine"               "Chloramphenicol" 
#>                             CTE                             CIP 
#>             "Chlortetracycline"                 "Ciprofloxacin" 
#>                             CLR                            CLA1 
#>                "Clarithromycin"               "Clavulanic acid" 
#>                             CLX                             CLI 
#>                 "Clinafloxacin"                   "Clindamycin" 
#>                             CLF                             CLO 
#>                   "Clofazimine"                   "Cloxacillin" 
#>                             COL                             CYC 
#>                      "Colistin"                   "Cycloserine" 
#>                             DAL                             DAN 
#>                   "Dalbavancin"                  "Danofloxacin" 
#>                             DAP                             DFX 
#>                    "Daptomycin"                  "Delafloxacin" 
#>                             DLM                             DIC 
#>                     "Delamanid"                 "Dicloxacillin" 
#>                             DIF                             DOR 
#>                    "Difloxacin"                     "Doripenem" 
#>                             DOX                             ENR 
#>                   "Doxycycline"                  "Enrofloxacin" 
#>                             ERV                             ETP 
#>                  "Eravacycline"                     "Ertapenem" 
#>                             ERY                             ETH 
#>                  "Erythromycin"                    "Ethambutol" 
#>                            ETI1                             FAR 
#>                   "Ethionamide"                     "Faropenem" 
#>                             FDX                             FLR 
#>                   "Fidaxomicin"                   "Florfenicol" 
#>                             FLC                             FLM 
#>                "Flucloxacillin"                    "Flumequine" 
#>                             FOS                             FUS 
#>                    "Fosfomycin"                  "Fusidic acid" 
#>                             GAM                             GAT 
#>                 "Gamithromycin"                  "Gatifloxacin" 
#>                             GEM                             GEN 
#>                  "Gemifloxacin"                    "Gentamicin" 
#>                             GEP                             IPM 
#>                   "Gepotidacin"                      "Imipenem" 
#>                             IMR                             INH 
#>           "Imipenem/relebactam"                     "Isoniazid" 
#>                             KAN                             LAS 
#>                     "Kanamycin"                     "Lasalocid" 
#>                             LMU                             LVX 
#>                     "Lefamulin"                  "Levofloxacin" 
#>                             LIN                             LNZ 
#>                    "Lincomycin"                     "Linezolid" 
#>                             MEC                             MEM 
#>                    "Mecillinam"                     "Meropenem" 
#>                             MEV                             MTR 
#>         "Meropenem/vaborbactam"                 "Metronidazole" 
#>                             MNO                             MON 
#>                   "Minocycline"               "Monensin sodium" 
#>                             MFX                             MUP 
#>                  "Moxifloxacin"                     "Mupirocin" 
#>                             NAL                             NAR 
#>                "Nalidixic acid"                       "Narasin" 
#>                             NEO                             NET 
#>                      "Neomycin"                    "Netilmicin" 
#>                             NIT                             NTR 
#>                "Nitrofurantoin"                   "Nitroxoline" 
#>                             NOR                             NVA 
#>                   "Norfloxacin"                 "Norvancomycin" 
#>                             OFX                             OMC 
#>                     "Ofloxacin"                  "Omadacycline" 
#>                             ORB                             ORI 
#>                  "Orbifloxacin"                   "Oritavancin" 
#>                             OXA                             OXO 
#>                     "Oxacillin"                 "Oxolinic acid" 
#>                             OXY                             PEF 
#>               "Oxytetracycline"                    "Pefloxacin" 
#>                             PHN                             PIP 
#>       "Phenoxymethylpenicillin"                  "Piperacillin" 
#>                             TZP                             PRL 
#>       "Piperacillin/tazobactam"                    "Pirlimycin" 
#>                             PRA                             PRI 
#>                 "Pradofloxacin"                 "Pristinamycin" 
#>                             PZA                             QDA 
#>                  "Pyrazinamide"     "Quinupristin/dalfopristin" 
#>                             RTP                             RZF 
#>                   "Retapamulin"                    "Rezafungin" 
#>                             RIB                             RIF 
#>                     "Rifabutin"                    "Rifampicin" 
#>                             RXT                             SAL 
#>                 "Roxithromycin"                   "Salinomycin" 
#>                             SEC                             SIT 
#>                   "Secnidazole"                  "Sitafloxacin" 
#>                             SPT                             SPI 
#>                 "Spectinomycin"                    "Spiramycin" 
#>                            STR1                             SUL 
#>                  "Streptomycin"                     "Sulbactam" 
#>                             SDI                             SMX 
#>                  "Sulfadiazine"              "Sulfamethoxazole" 
#>                             SOX                             TZD 
#>                 "Sulfisoxazole"                     "Tedizolid" 
#>                             TEC                             TLV 
#>                   "Teicoplanin"                    "Telavancin" 
#>                             TEM                             TCY 
#>                    "Temocillin"                  "Tetracycline" 
#>                             THI                             TIA 
#>                 "Thiamphenicol"                      "Tiamulin" 
#>                             TIC                             TCC 
#>                   "Ticarcillin"   "Ticarcillin/clavulanic acid" 
#>                             TGC                             TIP 
#>                   "Tigecycline"                  "Tildipirosin" 
#>                             TIL                             TOB 
#>                    "Tilmicosin"                    "Tobramycin" 
#>                             TMP                             SXT 
#>                  "Trimethoprim" "Trimethoprim/sulfamethoxazole" 
#>                             TUL                             TYL 
#>                 "Tulathromycin"                       "Tylosin" 
#>                            TYL1                             VAN 
#>                    "Tylvalosin"                    "Vancomycin" 
#>                             VIO 
#>                      "Viomycin" 
```
