#' Initial values for glmCensRd()
#'
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but other options are \code{"binomial"}, \code{"log-normal"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#'
#' @return A vector of parameters
#'
#' @export
#'
init_vals = function(Z = NULL, distY = "normal", distX = "normal") {
  if (distY %in% c("normal", "log-normal")) {
    params0 = c(rep(0, 2 + length(Z)), 1)
  } else if (distY == "binomial") {
    params0 = c(rep(0, 2 + length(Z)))
  } else if (distY %in% c("gamma", "inverse-gaussian", "weibull")) {
    params0 = rep(1E-4, 3 + length(Z))
  } else if (distY %in% c("exponential", "poisson")) {
    params0 = c(1E-4, rep(0, 1 + length(Z)))
  }

  if (distX %in% c("normal", "log-normal")) {
    params0 = c(params0, rep(0, 1 + length(Z)), 1)
  } else if (distX %in% c("gamma", "inverse-gaussian", "weibull")) {
    params0 = c(params0, rep(1E-4, 2 + length(Z)))
  } else if (distX %in% c("exponential", "poisson")) {
    params0 = c(params0, 1E-4, rep(0, length(Z)))
  }

  return(params0)
}
