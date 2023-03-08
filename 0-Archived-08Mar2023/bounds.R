#' Parameter bounds for glmCensRd()
#'
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but other options are \code{"binomial"}, \code{"log-normal"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#'
#' @return A list with the following elements:
#' \item{lower}{a vector of the lower bounds for all parameters.}
#' \item{upper}{a vector of the upper bounds for all parameters.}
#'
#' @export
#'
bounds = function(Z = NULL, distY = "normal", distX = "normal") {
  if (distY %in% c("normal", "log-normal")) {
    lower = c(rep(-Inf, 2 + length(Z)), 1E-4)
    upper = rep(Inf, 3 + length(Z))
  } else if (distY == "binomial") {
    lower = rep(-Inf, 2 + length(Z))
    upper = rep(Inf, 2 + length(Z))
  } else if (distY %in% c("gamma", "inverse-gaussian", "weibull")) {
    lower = rep(1E-4, 3 + length(Z))
    upper = rep(Inf, 3 + length(Z))
  } else if (distY %in% c("exponential", "poisson")) {
    lower = c(1E-4, rep(0, 1 + length(Z)))
    upper = rep(Inf, 2 + length(Z))
  }

  if (distX %in% c("normal", "log-normal")) {
    lower = c(lower, c(rep(-Inf, 1 + length(Z)), 0))
    upper = c(upper, rep(Inf, 2 + length(Z)))
  } else if (distX %in% c("gamma", "inverse-gaussian", "weibull")) {
    lower = c(lower, rep(1E-4, 2 + length(Z)))
    upper = c(upper, rep(Inf, 2 + length(Z)))
  } else if (distX %in% c("exponential", "poisson")) {
    params0 = c(params0, 1E-4, rep(0, length(Z)))
  }

  return(list(lower = lower, 
              upper = upper))
}
