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

```

## Logistic Regression
