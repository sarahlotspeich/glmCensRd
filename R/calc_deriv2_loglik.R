#' Calculate the second derivatives of the individual observations' log-likelihood contributions
#'
#' @param params Parameter values.
#' @param Y Name of outcome variable.
#' @param X Name of censored predictor variable.
#' @param W Name of observed (i.e., censored) version of \code{X}.
#' @param D Name of event indicator, defined to be = 1 if \code{X} was uncensored.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{NULL}
#' @param partX Size of partition of unobserved \code{X} for censored subjects. Default is \code{50}.
#' @param distY Distribution assumed for \code{Y} given \code{X} and \code{Z}.
#' @param distX Distribution assumed for \code{X} given \code{Z}.
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{W}, \code{D}, and \code{Z} (if applicable).
#' @export
#' @return A vector containing the derivatives of the log-likelihood contributions of the rows/observations in \code{data}.
#'
calc_deriv2_loglik <- function(params, Y, X, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data) {
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
    d_l0 <- calc_deriv_loglik(params = params0, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distY, distX = distX, data = data)

    # params + eps
    params1 <- matrix(data = params + eps * ej, nrow = p, ncol = 1)
    d_l1 <- calc_deriv_loglik(params = params1, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distY, distX = distX, data = data)

    # Estimate deriv = (l'(beta + eps) - l'(beta - eps)) / (2 * eps)
    d_theta[, start_col:(start_col + (p - 1))] <- (d_l1 - d_l0) / (2 * eps[j])
    start_col <- start_col + p
  }

  return(d_theta)
}
