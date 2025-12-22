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
#>                             AMC                             AMP 
#>   "Amoxicillin/clavulanic acid"                    "Ampicillin" 
#>                             SAM                             SAM 
#>          "Ampicillin/sulbactam"          "Ampicillin/sulbactam" 
#>                             APR                             APX 
#>                     "Apramycin"                  "Aspoxicillin" 
#>                             AVI                             AZM 
#>                    "Avilamycin"                  "Azithromycin" 
#>                             ATM                             AZA 
#>                     "Aztreonam"           "Aztreonam/avibactam" 
#>                             BAC                             BDQ 
#>                    "Bacitracin"                   "Bedaquiline" 
#>                             PEN                             CAP 
#>              "Benzylpenicillin"                   "Capreomycin" 
#>                             CEC                             CFR 
#>                      "Cefaclor"                    "Cefadroxil" 
#>                             LEX                             RID 
#>                     "Cefalexin"                  "Cefaloridine" 
#>                             CEP                             HAP 
#>                     "Cefalotin"                     "Cefapirin" 
#>                             CZO                             CDR 
#>                     "Cefazolin"                      "Cefdinir" 
#>                             FEP                             FPE 
#>                      "Cefepime"       "Cefepime/enmetazobactam" 
#>                             FPT                             FPZ 
#>           "Cefepime/tazobactam"           "Cefepime/zidebactam" 
#>                             FDC                             CFM 
#>                   "Cefiderocol"                      "Cefixime" 
#>                             CFP                             CSL 
#>                  "Cefoperazone"        "Cefoperazone/sulbactam" 
#>                             CSE                             CTX 
#>                     "Cefoselis"                    "Cefotaxime" 
#>                             CTT                             FOV 
#>                     "Cefotetan"                     "Cefovecin" 
#>                             FOX                             CPO 
#>                     "Cefoxitin"                     "Cefpirome" 
#>                             CPD                             CDC 
#>                   "Cefpodoxime"   "Cefpodoxime/clavulanic acid" 
#>                             CEQ                             CPT 
#>                    "Cefquinome"                   "Ceftaroline" 
#>                             CAZ                             CZA 
#>                   "Ceftazidime"         "Ceftazidime/avibactam" 
#>                             CTB                             TIO 
#>                    "Ceftibuten"                     "Ceftiofur" 
#>                             BPR                             CZT 
#>                  "Ceftobiprole"        "Ceftolozane/tazobactam" 
#>                             CZT                             CRO 
#>        "Ceftolozane/tazobactam"                   "Ceftriaxone" 
#>                             CXM                             CED 
#>                    "Cefuroxime"                    "Cephradine" 
#>                             CHL                             CTE 
#>               "Chloramphenicol"             "Chlortetracycline" 
#>                             CIP                             CLR 
#>                 "Ciprofloxacin"                "Clarithromycin" 
#>                            CLA1                             CLX 
#>               "Clavulanic acid"                 "Clinafloxacin" 
#>                             CLI                             CLF 
#>                   "Clindamycin"                   "Clofazimine" 
#>                             CLO                             COL 
#>                   "Cloxacillin"                      "Colistin" 
#>                             CYC                             DAL 
#>                   "Cycloserine"                   "Dalbavancin" 
#>                             DAN                             DAP 
#>                  "Danofloxacin"                    "Daptomycin" 
#>                             DFX                             DLM 
#>                  "Delafloxacin"                     "Delamanid" 
#>                             DIC                             DIF 
#>                 "Dicloxacillin"                    "Difloxacin" 
#>                             DOR                             DOX 
#>                     "Doripenem"                   "Doxycycline" 
#>                             ENR                             ERV 
#>                  "Enrofloxacin"                  "Eravacycline" 
#>                             ETP                             ERY 
#>                     "Ertapenem"                  "Erythromycin" 
#>                             ETH                            ETI1 
#>                    "Ethambutol"                   "Ethionamide" 
#>                             FAR                             FDX 
#>                     "Faropenem"                   "Fidaxomicin" 
#>                             FLR                             FLC 
#>                   "Florfenicol"                "Flucloxacillin" 
#>                             FLM                             FOS 
#>                    "Flumequine"                    "Fosfomycin" 
#>                             FUS                             GAM 
#>                  "Fusidic acid"                 "Gamithromycin" 
#>                             GAT                             GEM 
#>                  "Gatifloxacin"                  "Gemifloxacin" 
#>                             GEN                             GEP 
#>                    "Gentamicin"                   "Gepotidacin" 
#>                             IPM                             IMR 
#>                      "Imipenem"           "Imipenem/relebactam" 
#>                             INH                             KAN 
#>                     "Isoniazid"                     "Kanamycin" 
#>                             LAS                             LMU 
#>                     "Lasalocid"                     "Lefamulin" 
#>                             LVX                             LIN 
#>                  "Levofloxacin"                    "Lincomycin" 
#>                             LNZ                             MEC 
#>                     "Linezolid"                    "Mecillinam" 
#>                             MEM                             MEV 
#>                     "Meropenem"         "Meropenem/vaborbactam" 
#>                             MTR                             MNO 
#>                 "Metronidazole"                   "Minocycline" 
#>                             MON                             MFX 
#>               "Monensin sodium"                  "Moxifloxacin" 
#>                             MUP                             NAL 
#>                     "Mupirocin"                "Nalidixic acid" 
#>                             NAR                             NEO 
#>                       "Narasin"                      "Neomycin" 
#>                             NET                             NIT 
#>                    "Netilmicin"                "Nitrofurantoin" 
#>                             NTR                             NOR 
#>                   "Nitroxoline"                   "Norfloxacin" 
#>                             NVA                             OFX 
#>                 "Norvancomycin"                     "Ofloxacin" 
#>                             OMC                             ORB 
#>                  "Omadacycline"                  "Orbifloxacin" 
#>                             ORI                             OXA 
#>                   "Oritavancin"                     "Oxacillin" 
#>                             OXO                             OXY 
#>                 "Oxolinic acid"               "Oxytetracycline" 
#>                             PEF                             PHN 
#>                    "Pefloxacin"       "Phenoxymethylpenicillin" 
#>                             PIP                             TZP 
#>                  "Piperacillin"       "Piperacillin/tazobactam" 
#>                             PRL                             PRA 
#>                    "Pirlimycin"                 "Pradofloxacin" 
#>                             PRI                             PZA 
#>                 "Pristinamycin"                  "Pyrazinamide" 
#>                             QDA                             RTP 
#>     "Quinupristin/dalfopristin"                   "Retapamulin" 
#>                             RZF                             RIB 
#>                    "Rezafungin"                     "Rifabutin" 
#>                             RIF                             RXT 
#>                    "Rifampicin"                 "Roxithromycin" 
#>                             SAL                             SEC 
#>                   "Salinomycin"                   "Secnidazole" 
#>                             SIT                             SPT 
#>                  "Sitafloxacin"                 "Spectinomycin" 
#>                             SPI                            STR1 
#>                    "Spiramycin"                  "Streptomycin" 
#>                             SUL                             SDI 
#>                     "Sulbactam"                  "Sulfadiazine" 
#>                             SMX                             SOX 
#>              "Sulfamethoxazole"                 "Sulfisoxazole" 
#>                             TZD                             TEC 
#>                     "Tedizolid"                   "Teicoplanin" 
#>                             TLV                             TEM 
#>                    "Telavancin"                    "Temocillin" 
#>                             TCY                             THI 
#>                  "Tetracycline"                 "Thiamphenicol" 
#>                             TIA                             TIC 
#>                      "Tiamulin"                   "Ticarcillin" 
#>                             TCC                             TGC 
#>   "Ticarcillin/clavulanic acid"                   "Tigecycline" 
#>                             TIP                             TIL 
#>                  "Tildipirosin"                    "Tilmicosin" 
#>                             TOB                             TMP 
#>                    "Tobramycin"                  "Trimethoprim" 
#>                             SXT                             TUL 
#> "Trimethoprim/sulfamethoxazole"                 "Tulathromycin" 
#>                             TYL                            TYL1 
#>                       "Tylosin"                    "Tylvalosin" 
#>                             VAN                             VIO 
#>                    "Vancomycin"                      "Viomycin" 
```
