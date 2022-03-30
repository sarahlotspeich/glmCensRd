#' Calculate P(X|Z)
#'
#' @param x Predictor values (scalar or vector).
#' @param z Covariate values (scalar, vector, or dataframe).
#' @param distX Distribution assumed for \code{x} given \code{z}.
#' @param eta_params Vector of model parameters.
#'
#' @return A scalar or numeric vector the same length as the data input
#'
#' @export
#'
#'
calc_pXgivZ <- function(x, z = NULL, distX, eta_params) {
  if (distX %in% c("normal", "log-normal")) {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanX <- eta_params[1]
    if (!is.null(z)) {
      eta1 <- eta_params[-c(1, length(eta_params))]
      if (length(eta1) == 1) {
        meanX <- meanX + eta1 * z
      } else {
        meanX <- meanX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
      }
    }
    ## Estimate sqrt(variance) directly --------------
    sigX <- eta_params[length(eta_params)]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    eX <- x - meanX
    if (distX == "normal") {
      pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    } else {
      pXgivZ <- (1 / (x * sigX * sqrt(2 * pi))) * exp(- (log(x) - meanX) ^ 2 / (2 * sigX ^ 2))
    }
    # -------------------------------------- Calculate
  } else if (distX == "gamma") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- eta_params[1]
    ## Construct mean --------------------------------
    meanX <- eta_params[2]
    if (!is.null(z)) {
      eta1 <- eta_params[-c(1:2)]
      if (length(eta1) == 1) {
        meanX <- meanX + eta1 * z
      } else {
        meanX <- meanX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
      }
    }
    ## Construct scale -------------------------------
    scaleX <- meanX / shapeX
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- (1 / gamma(shapeX)) * scaleX ^ (- shapeX) * (x ^ (shapeX - 1)) * exp(- x / scaleX)
    # -------------------------------------- Calculate
  } else if (distX == "inverse-gaussian") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- eta_params[1]
    ## Construct the mean ----------------------------
    meanX <- eta_params[2]
    if (!is.null(z)) {
      eta1 <- eta_params[-c(1:2)]
      if (length(eta1) == 1) {
        meanX <- meanX + eta1 * z
      } else {
        meanX <- meanX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
      }
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- sqrt((shapeX / (2 * pi * x ^ 3))) * exp(- (shapeX * (x - meanX) ^ 2) / (2 * meanX ^ 2 * x))
    # -------------------------------------- Calculate
  } else if (distX == "weibull") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- eta_params[1]
    ## Construct scale -------------------------------
    scaleX <- eta_params[2]
    if (!is.null(z)) {
      eta1 <- eta_params[-c(1:2)]
      if (length(eta1) == 1) {
        scaleX <- scaleX + eta1 * z
      } else {
        scaleX <- scaleX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
      }
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- (shapeX / scaleX) * ((x / scaleX) ^ (shapeX - 1)) * exp(- (x / scaleX) ^ shapeX)
    # -------------------------------------- Calculate
  } else if (distX %in% c("exponential", "poisson")) {
    # Get parameters ---------------------------------
    ## Construct rate  -------------------------------
    rateX <- eta_params[1]
    if (!is.null(z)) {
      eta1 <- eta_params[-c(1)]
      if (length(eta1) == 1) {
        rateX <- rateX + eta1 * z
      } else {
        rateX <- rateX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
      }
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    if (distX == "exponential") {
      pXgivZ <- rateX * exp(- rateX * x)
    } else {
      pXgivZ <- rateX ^ x * exp(- rateX) / factorial(x)
    }
    # -------------------------------------- Calculate
  }
  return(pXgivZ)
}
