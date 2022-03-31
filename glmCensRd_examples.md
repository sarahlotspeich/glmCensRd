glmCensRd_examples
================
Sarah Lotspeich
3/3/2022

``` r
source("~sarahlotspeich/Dropbox/UNC/Sieve-Likelihood-Research/generate_data.R")
```

``` r
# Reproducibility
set.seed(333)
```

# `distY = "normal"`

## Linear Regression Outcome Model

### `distX = "normal"`

#### Predictor follows a normal distribution

``` r
# Generate data
temp <- generate_data(n = 1000, censoring = "light", distY = "normal", distX = "normal")
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub> = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *μ*<sub>*X*</sub> = 0.75 + 0.25*Z*, and
-   *σ*<sub>*X*</sub><sup>2</sup> = 1.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "normal", data = temp, steptol = 1e-2, iterlim = 100)
```

``` r
# Check convergence status
with(fit, code)
```

    ## [1] 2

The `code` is returned from `glmCensRd()`’s internal call to `nlm()` and
per the documentation is interpreted as follows:

-   `1`: relative gradient is close to zero, current iterate is probably
    solution.
-   `2`: successive iterates within tolerance, current iterate is
    probably solution.
-   `3`: last global step failed to locate a point lower than estimate.
    Either estimate is an approximate local minimum of the function or
    steptol is too small.
-   `4`: iteration limit exceeded.
-   `5`: maximum step size stepmax exceeded five consecutive times.
    Either the function is unbounded below, becomes asymptotic to a
    finite value from above in some direction or stepmax is too small.

In general, we accept `code <= 2` as successful convergence to the MLE.
Results of `code = 3-4` might require tuning the `steptol` or `iterlim`
parameters in the call to `glmCensRd()`.

``` r
# Inspect output 
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

    ## $distY
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.9961464 NA
    ## X           0.5024525 NA
    ## z           0.2635431 NA
    ## 
    ## $sigma2
    ## [1] 1.052978

From the `outcome_model` slot of the model `fit`, we get the following:

-   *μ̂*<sub>*Y*</sub>(*X*,*Z*) = 1.02 + 0.50*X* + 0.30*Z* and
-   *σ̂*<sub>*Y*</sub><sup>2</sup> = 1.05.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.7566012 NA
    ## z           0.2619449 NA
    ## 
    ## $sigma2
    ## [1] 1.003771

Then we get the remaining parameter estimates from `predictor_model`:

-   *μ̂*<sub>*X*</sub>(*Z*) = 0.76 + 0.31*Z* and
-   *σ̂*<sub>*X*</sub><sup>2</sup> = 0.98.

### `distX = "gamma"`

#### Predictor follows a gamma distribution

``` r
# Generate data
temp <- generate_data(n = 1000, censoring = "light", distY = "normal", distX = "gamma")
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *μ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*,
-   *ν*<sub>*X*</sub> = 2, and
-   *λ*<sub>*X*</sub>(*Z*) = (0.75+0.25*Z*)/2,

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "gamma", data = temp)
```

``` r
# Check convergence status
with(fit, code)
```

    ## [1] 2

``` r
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

    ## $distY
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 1.0512975 NA
    ## X           0.4375592 NA
    ## z           0.2591808 NA
    ## 
    ## $sigma2
    ## [1] 1.026661

The `code` and `outcome_model` slots are interpreted in the same way as
the previous section, since once again `distY = "normal"`.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "gamma"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.7264838 NA
    ## z           0.2361227 NA
    ## 
    ## $shape
    ##                coeff se
    ## (Intercept) 2.000851 NA

Then we get the remaining parameter estimates from `predictor_model`:

-   *μ̂*<sub>*X*</sub> = 0.78 + 0.24*Z*,
-   $\\hat{\\nu\_{X}} = 2.07$, and
-   $\\hat{\\lambda\_{X}} = \\hat{\\mu}\_X / \\hat{\\nu\_{X}} = \\left(0.78 + 0.24 Z\\right) / 2.07$.

### `distX = "inverse-gaussian"`

#### Predictor follows an inverse-Gaussian distribution

``` r
# Generate data
temp <- generate_data(n = 1000, censoring = "light", distY = "normal", distX = "inverse-gaussian")
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *μ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*, and
-   *λ*<sub>*X*</sub> = 2.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "inverse-gaussian", data = temp)
```

``` r
# Check convergence status
with(fit, code)
```

    ## [1] 2

``` r
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

    ## $distY
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.9868984 NA
    ## X           0.4845857 NA
    ## z           0.2764057 NA
    ## 
    ## $sigma2
    ## [1] 1.009602

The `code` and `outcome_model` slots are interpreted in the same way as
the previous section, since once again `distY = "normal"`.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "inverse-gaussian"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.7666817 NA
    ## z           0.2415719 NA
    ## 
    ## $shape
    ##                coeff se
    ## (Intercept) 2.069311 NA

Then we get the remaining parameter estimates from `predictor_model`:

-   *μ̂*<sub>*X*</sub> = 0.77 + 0.23*Z* and
-   $\\hat{\\lambda\_{X}} = 2.04$.

### `distX = "weibull"`

#### Predictor follows a Weibull distribution

``` r
# Generate data
temp <- generate_data(n = 1000, censoring = "light", distY = "normal", distX = "weibull")
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *λ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*, and
-   *ν*<sub>*X*</sub> = 2.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "weibull", data = temp)
```

``` r
# Check convergence status
with(fit, code)
```

    ## [1] 2

``` r
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

    ## $distY
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.9333027 NA
    ## X           0.5471663 NA
    ## z           0.2665987 NA
    ## 
    ## $sigma2
    ## [1] 1.006658

The `code` and `outcome_model` slots are interpreted in the same way as
the previous section, since once again `distY = "normal"`.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "weibull"
    ## 
    ## $scale
    ##                 coeff se
    ## (Intercept) 0.7974696 NA
    ## z           0.2289564 NA
    ## 
    ## $shape
    ##                coeff se
    ## (Intercept) 1.934245 NA

Then we get the remaining parameter estimates from `predictor_model`:

-   $\\hat{\\lambda\_{X}}\_X = 0.71 + 0.26 Z$ and
-   $\\hat{\\nu\_{X}} = 1.96$.

### `distX = "exponential"`

#### Predictor follows a exponential distribution

``` r
# Generate data
temp <- generate_data(n = 1000, censoring = "light", distY = "normal", distX = "exponential")
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *γ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "exponential", data = temp)
```

``` r
# Check convergence status
with(fit, code)
```

    ## [1] 2

``` r
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

    ## $distY
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.9559316 NA
    ## X           0.5128030 NA
    ## z           0.2721641 NA
    ## 
    ## $sigma2
    ## [1] 1.003454

The `code` and `outcome_model` slots are interpreted in the same way as
the previous section, since once again `distY = "normal"`.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "exponential"
    ## 
    ## $rate
    ##                 coeff se
    ## (Intercept) 0.7236389 NA
    ## z           0.2976122 NA

Then we get the remaining parameter estimates from `predictor_model`:

-   $\\hat{\\gamma\_{X}}(Z) = 0.74 + 0.28 Z$.

### `distX = "poisson"`

#### Predictor follows a Poisson distribution

``` r
# Generate data
temp <- generate_data(n = 1000, censoring = "light", distY = "normal", distX = "poisson")
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *γ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "poisson", data = temp)
```

``` r
# Check convergence status
with(fit, code)
```

    ## [1] 2

``` r
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

    ## $distY
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff se
    ## (Intercept) 0.8427277 NA
    ## X           0.5151241 NA
    ## z           0.2815058 NA
    ## 
    ## $sigma2
    ## [1] 0.9712626

The `code` and `outcome_model` slots are interpreted in the same way as
the previous section, since once again `distY = "normal"`.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "poisson"
    ## 
    ## $rate
    ##                 coeff se
    ## (Intercept) 0.7282190 NA
    ## z           0.1959749 NA

Then we get the remaining parameter estimates from `predictor_model`:

-   $\\hat{\\gamma\_{X}}(Z) = 0.64 + 0.43 Z$.
