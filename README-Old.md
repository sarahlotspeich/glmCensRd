# `glmCensRd`: Fitting generalized linear models with censored predictors in `R`

## Installation

```{r}
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/glmCensRd", ref = "main")
```

# Simulated Examples 
## Normal Linear Regression Outcome Model (`distY = "normal"`)

### `distX = "normal"`
#### Predictor follows a conditional normal distribution

```{r}
# Generate data
set.seed(333)
n <- 1000
z <- rlnorm(n = n)
x <- rnorm(n = n, mean = 0.75 + 0.25 * z, sd = 1)
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
c <- rweibull(n = n, shape = 1, scale = 3)
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
temp <- cbind(y, w, x, z, d)
```

This data generation scheme leads to true values: 

\begin{align*}
P_{\pmb{\theta}}(Y,X|Z) &= P_{\pmb{\beta}}(Y|X,Z)P_{\pmb{\eta}}(X|Z) \\
&= \frac{1}{\sqrt{2\pi\sigma_{Y}^2}}\exp\left\{-\frac{\left(Y - \mu_{Y}(X,Z)\right)^2}{2\sigma_{Y}^2}\right\}\frac{1}{\sqrt{2\pi\sigma_{X}^2}}\exp\left\{-\frac{\left(X - \mu_{X}(Z)\right)^2}{2\sigma_{X}^2}\right\}
\end{align*}
where 

  -   $\mu_Y = 1 + 0.5 X + 0.25 Z$, 
  -   $\sigma_{Y}^2 =1$, 
  -   $\mu_X = 0.75 + 0.25 Z$, and 
  -   $\sigma_{X}^2 =1$. 

```{r}
# Fit model 
fit <- glmCensRd::glmCensRd(Y = "y", W = "w", D = "d", Z = "z", distY = "normal", distX = "normal", data = temp)
```

```{r}
# Check convergence status
with(fit, code)
```

The `code` is returned from `glmCensRd()`'s internal call to `nlm()` and per the documentation is interpreted as follows:

  -   `1`: relative gradient is close to zero, current iterate is probably solution.
  -   `2`: successive iterates within tolerance, current iterate is probably solution.
  -   `3`: last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.
  -   `4`: iteration limit exceeded.
  -   `5`: maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.
  
In general, we accept `code <= 2` as successful convergence to the MLE. Results of `code = 3-4` might require tuning the `steptol` or `iterlim` parameters in the call to `glmCensRd()`. 

```{r}
# Inspect output 
## Outcome model P(Y|X,Z)
with(fit, outcome_model)
```

From the `outcome_model` slot of the model `fit`, we get the following:

  -   $\hat{\mu}_Y(X,Z) = 1.02 + 0.50 X + 0.30 Z$ and
  -   $\hat{\sigma}_{Y}^2 =1.05$.

```{r}
## Predictor model P(X|Z)
with(fit, predictor_model)
```

Then we get the remaining parameter estimates from `predictor_model`: 

  -   $\hat{\mu}_X(Z) = 0.76 + 0.31 Z$ and 
  -   $\hat{\sigma}_{X}^2 = 0.98$.

```{r}
# Fit model 
init_params <- c(0, 0, 0, var(y), 0, 0, var(x, na.rm = TRUE))
fit <- glmCensRd::glmCensRd(params0 = init_params, 
                            Y = "y", X = "x", W = "w", D = "d", Z = "z", 
                            distY = "normal", distX = "normal", data = temp)

# Inspect output
## Analysis model P(Y|X,Z)
fit$modY
$distY
[1] "normal"

$coeff
                coeff se
(Intercept) 1.0084869 NA
x           0.5117100 NA
z           0.2507936 NA

$sigma2
[1] 1.042889

## Predictor model P(X|Z)
fit$modX
$distX
[1] "normal"

$coeff
                coeff se
(Intercept) 0.7435398 NA
z           0.2004321 NA

$shape
[1] NA

$sigma2
[1] 1.046869

## Convergence code
fit$code
[1] 2
```

## Logistic Regression
