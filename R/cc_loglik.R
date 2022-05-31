#' Complete-case log-likelihood

#' @param params parameter values.
#' @param Y name of outcome variable.
#' @param X name of censored predictor variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data a dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but other options are \code{"binomial"}, \code{"log-normal"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"inverse-gaussian"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#'
#' @export
#'
#' @return A scalar of the log-likelihood function (negated for use with \code{nlm}()).
#'
cc_loglik <- function(params, Y, X, W, Z = NULL, data, distY = "normal", distX = "normal") {
  ####################################################
  # Joint density P(Y,X,Z) ###########################
  ####################################################
  pYXandZ <- calc_pYXandZ(y = data[, Y],
                          x = data[, X],
                          z = data[, Z],
                          distY = distY,
                          distX = distX,
                          params = params)

  # If params are out of domain, calc_pYXandZ returns 1E8
  if (any(is.na(pYXandZ))) {
    return(1E8)
  }

  ####################################################
  # Likelihood #######################################
  ####################################################
  ll <- sum(log(pYXandZ))

  # Return (-1) x log-likelihood for use with nlm() --
  return(- ll)
}
