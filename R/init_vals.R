init_vals <- function(Z, distY, distX) {
  if (distY %in% c("normal", "log-normal")) {
    params0 <- c(rep(0, length(c(1, 1, Z))), 1)
  } else if (distY == "binomial") {
    params0 <- c(rep(0, length(c(1, 1, Z))))
  }

  if (distX %in% c("normal", "log-normal")) {
    params0 <- c(params0, rep(0, length(c(1, Z))), 1)
  } else if (distX %in% c('gamma', "inverse-gaussian")) {
    params0 <- c(params0, 1E-4, rep(1E-4, length(c(1, Z))))
  } else if (distX == "weibull") {
    params0 <- c(params0, 1E-4, rep(1E-4, length(c(1, Z))))
  } else if (distX %in% c("exponential", "poisson")) {
    params0 <- c(params0, rep(1E-4, length(c(1, Z))))
  }

  return(params0)
}
