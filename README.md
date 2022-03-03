# `glmCensRd`: Fitting generalized linear models with censored predictors in `R`

## Installation

```{r}
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/glmCensRd", ref = "main")
```

# Simulated Examples 

## Normal Linear Regression 

### `distX = "normal"`

```{r}
# Generate data
n <- 1000
z <- rnorm(n = n)
x <- rnorm(n = n, mean = 0.75 + 0.25 * z, sd = 1)
y <- rnorm(n = n, mean = 1 + 0.5 * x + 0.25 * z, sd = 1)
c <- rweibull(n = n, shape = 1, scale = 3)
w <- pmin(x, c)
d <- as.numeric(x <= c)
x[d == 0] <- NA
temp <- cbind(y, w, x, z, d)

# Fit model 
init_params <- c(0, 0, 0, var(y), 0, 0, var(x, na.rm = TRUE))
fit <- glmCensRd::glmCensRd(params0 = init_params, 
                            Y = "y", X = "x", W = "w", D = "d", Z = "z", 
                            distY = "normal", distX = "normal", data = temp)

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

fit$code
[1] 2
```

## Logistic Regression
