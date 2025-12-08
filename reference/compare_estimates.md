# Plot to Compare Two Sets of Estimates

This function compares two sets of estimates by creating a plot that
overlays the estimates and confidence intervals for both sets. It can
also display the estimates in two separate plots.

## Usage

``` r
compare_estimates(
  tbl1,
  tbl2,
  title1 = NULL,
  title2 = NULL,
  title = NULL,
  sig = 0.05,
  colors = c("maroon", "blue4"),
  x_title = "Coefficient (95% CI)",
  y_title = "Variant",
  axis_label_size = 9,
  single_plot = TRUE,
  pd = position_dodge(width = 0.8),
  marker_order = NULL
)
```

## Arguments

- tbl1:

  A tibble containing the first set of summary statistics (e.g.,
  coefficients, p-values, CI) for each variant. Expected columns are:

  - `marker`: The name of the marker (e.g., variable name).

  - `pval`: The p-value for each marker.

  - `ci.lower`: The lower bound of the confidence interval.

  - `ci.upper`: The upper bound of the confidence interval.

  - `est`: The estimated coefficient.

- tbl2:

  A tibble containing the second set of summary statistics for each
  variant (same format as tbl1).

- title1:

  Title for tbl1 data. If single_plot, this will be the legend label for
  tbl1 data; otherwise it will be the title for the tbl1 plot.

- title2:

  Title for tbl2 data. If single_plot, this will be the legend label for
  tbl2 data; otherwise it will be the title for the tbl2 plot.

- title:

  (optional) The main title of the combined plot, if single_plot is
  TRUE.

- sig:

  (optional) For individual plots, the p-value threshold for data points
  to highlight as significant. Defaults to 0.05.

- colors:

  (optional) For combined plot, a vector of two colors to represent the
  two input datasets.

- x_title:

  (optional) The title for the x-axis. Defaults to "Coefficient (95%
  CI)".

- y_title:

  (optional) The title for the y-axis. Defaults to "Variant".

- axis_label_size:

  (optional) The font size of the axis labels. Defaults to 9.

- single_plot:

  A boolean indicating whether to make a single combined plot (TRUE), or
  plot each dataset side-by-side (FALSE).

- pd:

  (optional) Position dodge, i.e. spacing for the 2 estimates to be
  positioned above/below the line. Default 'position_dodge(width = 0.8)'

- marker_order:

  (optional) Vector indicating the order of the markers to be plotted on
  the y-axis.

## Value

A ggplot object displaying the comparison of the two sets of estimates.
