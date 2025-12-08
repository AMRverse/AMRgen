# Plot Combined Statistics of Logistic Regression and Solo PPV

This function creates a plot comparing logistic regression coefficients
with PPV values from solo marker analysis. It highlights markers based
on statistical significance of the logistic regression.

## Usage

``` r
plot_solo_logReg(combined_stats, sig = 0.05, title = NULL)
```

## Arguments

- combined_stats:

  A data frame containing combined statistics from logistic regression
  and solo PPV analysis.

- sig:

  A numeric value specifying the significance threshold for p-values.
  Default is 0.05.

- title:

  An optional character string specifying the plot title.

## Value

A ggplot2 object visualizing the relationship between PPV and logistic
regression estimates, with confidence intervals and significance
annotation.
