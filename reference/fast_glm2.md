# Run many regression models quickly, for many response variables.

Performs any type of cross-sectional regression analysis quickly for
many response variables.

## Usage

``` r
fast_glm2(
  response = NULL,
  response_name = NULL,
  Xvar = NULL,
  dataset = NULL,
  adj_vars = NULL,
  family = "gaussian",
  exponentiate = FALSE,
  p_adjust_meth = "fdr",
  robust = FALSE,
  sort = TRUE,
  decimals = 3
)
```

## Arguments

- response:

  Character vector of variable names to use as response of interest in
  separate regression models.

- response_name:

  Label to display for each response of interest

- Xvar:

  Main predictor to be used in each regression.

- dataset:

  Dataset in which these variables exist.

- adj_vars:

  Character vector of adjustment covariate names.

- family:

  Family for regression model, e.g. 'gaussian', 'binomial', 'poisson'.

- exponentiate:

  Whether or not to exponentiate the estimate from the model (e.g. for
  logistic or Poisson regression).

- p_adjust_meth:

  Method by which to provide adjustment of P-values.

- robust:

  If TRUE, returns robust standard error estimates.

- sort:

  If TRUE, sorts by p-value.

- decimals:

  Specifies the number of decimals to which to round estimates.
