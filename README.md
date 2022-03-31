`glmCensRd`: Fitting generalized linear models with censored predictors
in `R`
================
Sarah C. Lotspeich, Peter Guan, and Tanya P. Garcia

``` r
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/glmCensRd", ref = "main")
```

``` r
# Reproducibility
set.seed(333)
```

# Linear Regression Outcome Model

## Set `distY = "normal"`

### Predictor follows a normal distribution

#### Set `distX = "normal"`

``` r
# Generate data
n <- 1000
## Generate fully observed covariate
z <- rlnorm(n = n)
## Generate the (censored) predictor
x <- rnorm(n = n, mean = 0.75 + 0.25 * z, sd = 1)  
## Generate censoring value
c <- rweibull(n = n, shape = 1, scale = 3)  
## Generate the outcome 
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
## Calculate observed predictor + event indicator
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
# Construct dataframe 
temp <- data.frame(y, w, x, z, d)
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub> = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *μ*<sub>*X*</sub> = 0.75 + 0.25*Z*, and
-   *σ*<sub>*X*</sub><sup>2</sup> = 1.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "normal", data = temp, steptol = 1e-2, iterlim = 100, estSE = TRUE)
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
    ##                 coeff         se
    ## (Intercept) 0.9730323 0.09091957
    ## X           0.4566375 0.08317825
    ## z           0.2815724 0.02841242
    ## 
    ## $sigma2
    ## [1] 0.9893482

From the `outcome_model` slot of the model `fit`, we get the following:

-   *μ̂*<sub>*Y*</sub>(*X*,*Z*) = 0.973 + 0.457*X* + 0.282*Z* and
-   *σ̂*<sub>*Y*</sub><sup>2</sup> = 0.989.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff         se
    ## (Intercept) 0.7033991 0.05449020
    ## z           0.2993972 0.04264013
    ## 
    ## $sigma2
    ## [1] 1.000342

Then we get the remaining parameter estimates from `predictor_model`:

-   *μ̂*<sub>*X*</sub>(*Z*) = 0.703 + 0.299*Z* and
-   *σ̂*<sub>*X*</sub><sup>2</sup> = 1.000.

### Predictor follows a gamma distribution

#### Set `distX = "gamma"`

``` r
# Generate data
n <- 1000
## Generate fully observed covariate
z <- rlnorm(n = n)
## Generate the (censored) predictor
x <- rgamma(n = n, shape = 2, scale = (0.75 + 0.25 * z) / 2)
## Generate censoring value
c <- rweibull(n = n, shape = 1, scale = 3)  
## Generate the outcome 
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
## Calculate observed predictor + event indicator
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
# Construct dataframe 
temp <- data.frame(y, w, x, z, d)
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
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "gamma", data = temp, estSE = TRUE)
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
    ##                 coeff         se
    ## (Intercept) 1.0323830 0.04597472
    ## X           0.5267726 0.02059062
    ## z           0.2339656 0.01716312
    ## 
    ## $sigma2
    ## [1] 1.001325

The `code` and `outcome_model` slots are interpreted in the same way as
the first section, since once again `distY = "normal"`, with

-   *μ̂*<sub>*Y*</sub>(*X*,*Z*) = 1.032 + 0.527*X* + 0.234*Z* and
-   *σ̂*<sub>*Y*</sub><sup>2</sup> = 1.001.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "gamma"
    ## 
    ## $mean
    ##                 coeff         se
    ## (Intercept) 0.7075069 0.02304275
    ## z           0.2458071 0.02226052
    ## 
    ## $shape
    ##                coeff         se
    ## (Intercept) 2.027891 0.06915529

Then we get the remaining parameter estimates from `predictor_model`:

-   *μ̂*<sub>*X*</sub> = 0.708 + 0.246*Z*,
-   *ν̂*<sub>*X*</sub> = 2.028, and
-   *λ̂*<sub>*X*</sub> = *μ̂*<sub>*X*</sub>/*ν̂*<sub>*X*</sub> = (0.708+0.246*Z*)/2.028.

### Predictor follows an inverse-Gaussian distribution

#### Set `distX = "inverse-gaussian"`

``` r
# Generate data
n <- 1000
## Generate fully observed covariate
z <- rlnorm(n = n)
## Generate the (censored) predictor
x <- SuppDists::rinvGauss(n = n, nu = 0.75 + 0.25 * z, lambda = 2)
## Generate censoring value
c <- rweibull(n = n, shape = 1, scale = 3)  
## Generate the outcome 
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
## Calculate observed predictor + event indicator
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
# Construct dataframe 
temp <- data.frame(y, w, x, z, d)
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *μ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*, and
-   *λ*<sub>*X*</sub> = 2.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "inverse-gaussian", data = temp, estSE = TRUE)
```

``` r
# Check convergence status
with(fit, code)
```

    ## [1] 1

``` r
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

    ## $distY
    ## [1] "normal"
    ## 
    ## $mean
    ##                 coeff          se
    ## (Intercept) 0.9769124 0.040886988
    ## X           0.5201065 0.003835241
    ## z           0.2435986 0.016672972
    ## 
    ## $sigma2
    ## [1] 0.9942551

The `code` and `outcome_model` slots are interpreted in the same way as
the first section, since once again `distY = "normal"`, with

-   *μ̂*<sub>*Y*</sub>(*X*,*Z*) = 0.977 + 0.520*X* + 0.243*Z* and
-   *σ̂*<sub>*Y*</sub><sup>2</sup> = 0.994.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "inverse-gaussian"
    ## 
    ## $mean
    ##                 coeff         se
    ## (Intercept) 0.7854897 0.03956616
    ## z           0.2288820 0.03414016
    ## 
    ## $shape
    ##                coeff         se
    ## (Intercept) 2.031719 0.09624681

Then we get the remaining parameter estimates from `predictor_model`:

-   *μ̂*<sub>*X*</sub> = 0.785 + 0.229*Z* and
-   *λ̂*<sub>*X*</sub> = 2.032.

### Predictor follows a Weibull distribution

#### Set `distX = "weibull"`

``` r
# Generate data
n <- 1000
## Generate fully observed covariate
z <- rlnorm(n = n)
## Generate the (censored) predictor
x <- rweibull(n = n, shape = 2, scale = 0.75 + 0.25 * z)
## Generate censoring value
c <- rweibull(n = n, shape = 1, scale = 3)  
## Generate the outcome 
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
## Calculate observed predictor + event indicator
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
# Construct dataframe 
temp <- data.frame(y, w, x, z, d)
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *λ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*, and
-   *ν*<sub>*X*</sub> = 2.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "weibull", data = temp, estSE = TRUE)
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
    ##                 coeff         se
    ## (Intercept) 0.8642203 0.06081294
    ## X           0.5924190 0.05924627
    ## z           0.2713237 0.01815378
    ## 
    ## $sigma2
    ## [1] 0.9295751

The `code` and `outcome_model` slots are interpreted in the same way as
the first section, since once again `distY = "normal"`, with

-   *μ̂*<sub>*Y*</sub>(*X*,*Z*) = 0.864 + 0.592*X* + 0.271*Z* and
-   *σ̂*<sub>*Y*</sub><sup>2</sup> = 0.930.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "weibull"
    ## 
    ## $scale
    ##                 coeff         se
    ## (Intercept) 0.8015621 0.01755491
    ## z           0.2238740 0.01529947
    ## 
    ## $shape
    ##                coeff         se
    ## (Intercept) 1.925981 0.04179089

Then we get the remaining parameter estimates from `predictor_model`:

-   *λ̂*<sub>*X*<sub>*X*</sub></sub> = 0.802 + 0.224*Z* and
-   *ν̂*<sub>*X*</sub> = 1.923.

### Predictor follows a exponential distribution

#### Set `distX = "exponential"`

``` r
# Generate data
n <- 1000
## Generate fully observed covariate
z <- rlnorm(n = n)
## Generate the (censored) predictor
x <- rexp(n = n, rate = 0.75 + 0.25 * z)
## Generate censoring value
c <- rweibull(n = n, shape = 1, scale = 3)  
## Generate the outcome 
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
## Calculate observed predictor + event indicator
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
# Construct dataframe 
temp <- data.frame(y, w, x, z, d)
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *γ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "exponential", data = temp, estSE = TRUE)
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
    ##                 coeff         se
    ## (Intercept) 1.1186313 0.11496119
    ## X           0.4680319 0.01892211
    ## z           0.2151738 0.06292374
    ## 
    ## $sigma2
    ## [1] 0.9376589

The `code` and `outcome_model` slots are interpreted in the same way as
the first section, since once again `distY = "normal"`, with

-   *μ̂*<sub>*Y*</sub>(*X*,*Z*) = 1.119 + 0.468*X* + 0.215*Z* and
-   *σ̂*<sub>*Y*</sub><sup>2</sup> = 0.934.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "exponential"
    ## 
    ## $rate
    ##                 coeff         se
    ## (Intercept) 0.7817447 0.07201083
    ## z           0.2815443 0.06723741

Then we get the remaining parameter estimates from `predictor_model`:

-   *γ̂*<sub>*X*</sub>(*Z*) = 0.782 + 0.282*Z*.

### Predictor follows a Poisson distribution

#### Set `distX = "poisson"`

``` r
# Generate data
n <- 1000
## Generate fully observed covariate
z <- rlnorm(n = n)
## Generate the (censored) predictor
x <- rpois(n = n, lambda = 0.75 + 0.25 * z)
## Generate censoring value
c <- rweibull(n = n, shape = 1, scale = 3)  
## Generate the outcome 
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
## Calculate observed predictor + event indicator
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
# Construct dataframe 
temp <- data.frame(y, w, x, z, d)
```

This data generation scheme leads to true values:

where

-   *μ*<sub>*Y*</sub>(*X*,*Z*) = 1 + 0.5*X* + 0.25*Z*,
-   *σ*<sub>*Y*</sub><sup>2</sup> = 1,
-   *γ*<sub>*X*</sub>(*Z*) = 0.75 + 0.25*Z*.

``` r
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "poisson", data = temp, estSE = TRUE)
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
    ##                 coeff         se
    ## (Intercept) 0.9477904 0.05371840
    ## X           0.5084887 0.03455069
    ## z           0.2205852 0.02521094
    ## 
    ## $sigma2
    ## [1] 1.00549

The `code` and `outcome_model` slots are interpreted in the same way as
the first section, since once again `distY = "normal"`, with

-   *μ̂*<sub>*Y*</sub>(*X*,*Z*) = 0.948 + 0.508*X* + 0.221*Z* and
-   *σ̂*<sub>*Y*</sub><sup>2</sup> = 1.005.

``` r
## Predictor model P(X|Z)
with(fit, predictor_model)
```

    ## $distX
    ## [1] "poisson"
    ## 
    ## $rate
    ##                 coeff         se
    ## (Intercept) 0.7446170 0.05068968
    ## z           0.1911478 0.03251821

Then we get the remaining parameter estimates from `predictor_model`:

-   *γ̂*<sub>*X*</sub>(*Z*) = 0.745 + 0.191*Z*.
