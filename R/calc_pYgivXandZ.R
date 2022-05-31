#' Calculate P(Y|X,Z)
#'
#' @param y Outcome values (scalar or vector).
#' @param x Predictor values (scalar or vector).
#' @param z Covariate values (scalar, vector, or dataframe).
#' @param distY Distribution assumed for \code{y} given \code{x} and \code{z}.
#' @param beta_params Vector of model parameters.
#'
#' @return A scalar or numeric vector the same length as the data input
#'
#' @export
#'
#'
calc_pYgivXandZ <- function(y, x, z = NULL, distY, beta_params) {
  if (distY == "normal") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- beta_params[1] + beta_params[2] * matrix(data = x, ncol = 1)
    if (!is.null(z)) {
      beta2 <- beta_params[-c(1:2, length(beta_params))]
      if (length(beta2) == 1) {
        meanY <- meanY + beta2 * z
      } else {
        meanY <- meanY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
      }
    }
    ## Estimate sqrt(variance) directly --------------
    sigY <- beta_params[length(beta_params)]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    eY <- as.numeric(y) - meanY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # -------------------------------------- Calculate
  } else if (distY == "binomial") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- beta_params[1] + beta_params[2] * matrix(data = x, ncol = 1)
    if (!is.null(z)) {
      beta2 <- beta_params[-c(1:2)]
      if (length(beta2) == 1) {
        meanY <- meanY + z
      } else {
        meanY <- meanY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
      }
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    # pYgivXZ <- exp(- (1 - y) * meanY) / (1 + exp(meanY))
    # pYgivXZ <- dbinom(x = as.numeric(y), size = 1, prob = 1 / (1 + exp(- meanY)))
    # -------------------------------------- Calculate
  } else if (distY == "poisson") {
    # Get parameters -----------------------------------
    # ----------------------------------- Get parameters
    # Calculate ----------------------------------------
    muY <- exp(beta_params[1] + beta_params[2] * matrix(data = x, ncol = 1))
    if (!is.null(z)) {
      beta2 <- beta_params[-c(1:2)]
      if (length(beta2) == 1) {
        muY <- muY * exp(beta2 * z)
      } else {
      muY <- muY * exp(as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1)))
      }
    }
    pYgivXZ <- ((muY)^(as.numeric(y))) * exp(-muY) / factorial(as.numeric(y))
    # ---------------------------------------- Calculate
  }
  return(pYgivXZ)
}
