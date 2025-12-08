# Extract Details from a logistf Model

This function extracts and formats the estimates, confidence intervals,
and p-values from a fitted logistf model.

## Usage

``` r
logistf_details(model)
```

## Arguments

- model:

  A fitted logistf model object.

## Value

A tibble containing the estimates, confidence intervals, and p-values
for each predictor in the model.

## Examples

``` r
library(logistf)
library(ggplot2)
model <- logistf(case ~ age + oc + vic + vicl + vis + dia, data = sex2)
model_details <- logistf_details(model)
autoplot(model_details)
```
