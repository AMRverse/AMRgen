# Merge Logistic Regression and Solo PPV Statistics

This function merges logistic regression model statistics with solo PPV
statistics and creates a combined plot.

## Usage

``` r
merge_logreg_soloppv(model, solo_stats, title = NULL, plot = TRUE)
```

## Arguments

- model:

  A data frame containing logistic regression model statistics.

- solo_stats:

  A data frame containing solo PPV statistics.

- title:

  An optional title for the plot.

- plot:

  Logical indicating whether to generate the plot.

## Value

A list containing:

- `combined`: A merged data frame of logistic regression and PPV
  statistics.

- `plot`: A ggplot object showing the relationship between model
  estimates and PPV.

## Examples

``` r
if (FALSE) { # \dontrun{
soloPPV_cipro <- solo_ppv_analysis(ecoli_geno, ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno"
)
logistic_cipro <- amr_logistic(ecoli_geno, ecoli_ast,
  "Ciprofloxacin", c("Quinolones"),
  maf = 5
)
allstatsR <- merge_logreg_soloppv(logistic_cipro$modelR,
  soloPPV_cipro$solo_stats %>% filter(category == "R"),
  title = "Quinolone markers vs Cip R"
)
} # }
```
