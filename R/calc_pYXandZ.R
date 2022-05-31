#' Calculate P(Y,X,Z) = P(Y|X,Z)P(X|Z)
#'
#' @param y Outcome values (scalar or vector).
#' @param x Predictor values (scalar or vector).
#' @param z Covariate values (scalar, vector, or dataframe).
#' @param distY Distribution assumed for \code{y} given \code{x} and \code{z}.
#' @param beta_params Vector of model parameters.
#' @param distX Distribution assumed for \code{x} given \code{z}.
#' @param eta_params Vector of model parameters.
#'
#' @return A scalar or numeric vector the same length as the data input
#'
#' @export
#'
#'
calc_pYXandZ <- function(y, x, z = NULL, distY, beta_params, distX, eta_params) {
  # Calculate P(Y|X,Z)
  pYgivXandZ <- calc_pYgivXandZ(y = y,
                                x = x,
                                z = z,
                                distY = distY,
                                beta_params = beta_params)

  # Calculate P(X|Z)
  pXgivZ <- calc_pXgivZ(x = x,
                        z = z,
                        distX = distX,
                        eta_params = eta_params)

  # Return joint density
  return(pYgivXandZ * pXgivZ)
}
