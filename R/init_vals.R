init_vals <- function(Z, distY, distX) {
  if (distY %in% c("normal", "log-normal")) {
    params0 <- c(rep(0, 2 + length(Z)), 1)
  } else if (distY == "binomial") {
    params0 <- c(rep(0, 2 + length(Z)))
  } else if (distY %in% c("gamma", "inverse-gaussian", "weibull")) {
    params0 <- rep(1E-4, 3 + length(Z))
  } else if (distY %in% c("exponential", "poisson")) {
    params0 <- c(1E-4, rep(0, 1 + length(Z)))
  }

  if (distX %in% c("normal", "log-normal")) {
    params0 <- c(params0, rep(0, 1 + length(Z)), 1)
  } else if (distX %in% c("gamma", "inverse-gaussian", "weibull")) {
    params0 <- c(params0, 1E-4, rep(0, 1 + length(Z)))
  } else if (distX %in% c("exponential", "poisson")) {
    params0 <- c(params0, 1E-4, rep(0, length(Z)))
  }

  return(params0)
}
