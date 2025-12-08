# Extract Details from a Generalized Linear Model

This function extracts and formats the estimates, confidence intervals,
and p-values from a fitted glm model.

## Usage

``` r
glm_details(model)
```

## Arguments

- model:

  A fitted glm model object.

## Value

A tibble containing the estimates, confidence intervals, and p-values
for each predictor in the model.

## Examples

``` r
# Generate example data
set.seed(1)
dat <- data.frame(
  R = rbinom(100, 1, 0.3),
  geneA = rbinom(100, 1, 0.2),
  geneB = rbinom(100, 1, 0.5),
  geneC = rbinom(100, 1, 0.3)
)

# Fit logistic regression model and extract details
model <- glm(R ~ ., data = dat, family = binomial(link = "logit"))
model_details <- glm_details(model)
#> Waiting for profiling to be done...

# Plot model summary
ggplot2::autoplot(model_details)
```
