#' Calculate the second derivatives of the individual observations' log-likelihood contributions
#'
#' @param params parameter values.
#' @param Y name of outcome variable.
#' @param X name of censored predictor variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param D name of event indicator, defined to be \code{= 1} if \code{X} was uncensored and \code{0} otherwise.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data a dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param subdivisions (fed to \code{integrate}) the maximum number of subintervals used to integrate over unobserved \code{X} for censored subjects. Default is \code{100}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but \code{"binomial"} is the other option.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"inverse-gaussian"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @export
#' @return A vector containing the derivatives of the log-likelihood contributions of the rows/observations in \code{data}.
#'
calc_deriv2_loglik <- function(params, Y, X, W, D, Z = NULL, data, subdivisions = 100, distY = "normal", distX = "normal") {
  p <- length(params)
  eps <- params * (10 ^ (- 4))

  # Create matrix to save derivatives in
  d_theta <- matrix(data = 0, nrow = nrow(data), ncol = p ^ 2, byrow = FALSE)
  start_col <- 1
  for (j in 1:p) {
    # Create jth euclidean vector
    ej <- matrix(data = 0, nrow = p, ncol = 1)
    ej[j] <- 1

    # params - eps
    params0 <- matrix(data = params - eps * ej, nrow = p, ncol = 1)
    d_l0 <- calc_deriv_loglik(params = params0, Y = Y, X = X, W = W, D = D, Z = Z, subdivisions = subdivisions, distY = distY, distX = distX, data = data)

    # params + eps
    params1 <- matrix(data = params + eps * ej, nrow = p, ncol = 1)
    d_l1 <- calc_deriv_loglik(params = params1, Y = Y, X = X, W = W, D = D, Z = Z, subdivisions = subdivisions, distY = distY, distX = distX, data = data)

    # Estimate deriv = (l'(beta + eps) - l'(beta - eps)) / (2 * eps)
    d_theta[, start_col:(start_col + (p - 1))] <- (d_l1 - d_l0) / (2 * eps[j])
    start_col <- start_col + p
  }

  return(d_theta)
}
