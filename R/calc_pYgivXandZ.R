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
  if (distY %in% c("normal", "log-normal")) {
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
    sigY <- sqrt(beta_params[length(beta_params)] ^ 2)
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    # eY <- as.numeric(y) - meanY
    if (distY == "normal") {
      # pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
      pYgivXZ <- dnorm(x = y, mean = meanY, sd = sigY)
    } else {
      pYgivXZ <- dlnorm(x = y, meanlog = meanY, sdlog = sigY)
    }
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
    pYgivXZ <- dbinom(x = y, size = 1, prob = 1 / (1 + exp(- meanY)))
    pYgivXZ[y == 0] <- 1 - pYgivXZ[y == 0]
    # pYgivXZ <- exp(- (1 - y) * meanY) / (1 + exp(meanY))
    # pYgivXZ <- dbinom(x = as.numeric(y), size = 1, prob = 1 / (1 + exp(- meanY)))
    # -------------------------------------- Calculate
  } else if (distY == "weibull") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeY <- beta_params[1]
    ## Construct scale -------------------------------
    scaleY <- beta_params[2] + beta_params[3] * x
    if (!is.null(z)) {
      beta2 <- beta_params[-c(1:3)]
      if (length(eta1) == 1) {
        scaleY <- scaleY + beta2 * z
      } else {
        scaleY <- scaleY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
      }
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pYgivXZ <- dweibull(x = y, shape = shapeY, scale = scaleY)
    # -------------------------------------- Calculate
  } else if (distY %in% c("exponential", "poisson")) {
    # Get parameters ---------------------------------
    ## Construct rate  -------------------------------
    rateY <- beta_params[1] + beta_params[2] * x
    if (!is.null(z)) {
      beta2 <- beta_params[-c(1:2)]
      if (length(beta2) == 1) {
        rateY <- rateY + beta2 * z
      } else {
        rateY <- rateY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
      }
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    if (distX == "exponential") {
      pYgivXZ <- dexp(x = y, rate = rateY)
    } else {
      pYgivXZ <- dpois(x = y, lambda = rateY)
    }
    # -------------------------------------- Calculate
  }
  return(pYgivXZ)
}
