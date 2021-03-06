#' Complete-case log-likelihood

#' @param params parameter values.
#' @param Y name of outcome variable.
#' @param X name of censored predictor variable.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data a dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but other options are \code{"binomial"}, \code{"log-normal"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"inverse-gaussian"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#'
#' @export
#'
#' @return A scalar of the log-likelihood function (negated for use with \code{nlm}()).
#'
cc_loglik <- function(params, Y, X, Z = NULL, data, distY = "normal", distX = "normal") {
  ####################################################
  # Joint density P(Y,X,Z) ###########################
  ####################################################
  pYXandZ <- calc_pYXandZ(x = data[, X],
                          y = data[, Y],
                          z = data[, Z],
                          lengthZ = length(Z),
                          distY = distY,
                          distX = distX,
                          params = params)

  ####################################################
  # Likelihood #######################################
  ####################################################
  if (any(is.na(pYXandZ))) {
    # If params are out of domain, calc_pYXandZ returns NA
    ## And the log-likelihood needs to be arbitrarily "huge"
    return(1E8)
  } else {
    # Replace pYXandZ = 0 with 1 so that log(pYXandZ) = 0
    ## Otherwise log(pYXandZ) = -Inf
    pYXandZ[pYXandZ == 0] <- 1
    ll <- sum(log(pYXandZ))
    # Return (-1) x log-likelihood for use with nlm() --
    return(- ll)
  }
}
