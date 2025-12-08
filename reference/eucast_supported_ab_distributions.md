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
#>                             GEN                             IPM 
#>                    "Gentamicin"                      "Imipenem" 
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
