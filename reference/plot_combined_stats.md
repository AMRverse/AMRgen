# Plot Combined Statistics

This function creates a plot of combined logistic regression and solo
PPV statistics.

## Usage

``` r
plot_combined_stats(combined_stats, sig = 0.05, title = NULL)
```

## Arguments

- combined_stats:

  A data frame containing combined statistics from logistic regression
  and solo PPV.

- sig:

  A significance level for the logistic regression p-values. Default is
  0.05.

- title:

  An optional title for the plot.

## Value

A ggplot2 object representing the combined plot.
