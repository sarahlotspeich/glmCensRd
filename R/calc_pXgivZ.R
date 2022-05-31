#' Calculate P(X|Z)
#'
#' @param x Predictor values (scalar or vector).
#' @param z Covariate values (scalar, vector, or dataframe).
#' @param distX Distribution assumed for \code{x} given \code{z}.
#' @param eta_params Vector of model parameters.
#'
#' @importFrom SuppDists dinvGauss
#'
#' @export
#'
#' @return A scalar or numeric vector the same length as the data input
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
    sigX <- sqrt(eta_params[length(eta_params)] ^ 2)
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    # eX <- x - meanX
    if (distX == "normal") {
      # pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
      pXgivZ <- dnorm(x = x, mean = meanX, sd = sigX)
    } else {
      # pXgivZ <- (1 / (x * sigX * sqrt(2 * pi))) * exp(- (log(x) - meanX) ^ 2 / (2 * sigX ^ 2))
      pXgivZ <- dlnorm(x = x, meanlog = meanX, sdlog = sigX)
    }
    # -------------------------------------- Calculate
  } else if (distX == "gamma") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- exp(eta_params[1])
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
    # pXgivZ <- (1 / gamma(shapeX)) * scaleX ^ (- shapeX) * (x ^ (shapeX - 1)) * exp(- x / scaleX)
    suppressWarnings(
      pXgivZ <- dgamma(x = x, shape = shapeX, scale = scaleX)
    )
    # -------------------------------------- Calculate
    # Check: shape and scale of Gamma > 0 ------------
    pXgivZ[scaleX <= 0] <- NA
    # ------------ Check: shape and scale of Gamma > 0
  } else if (distX == "inverse-gaussian") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- exp(eta_params[1])
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
    # pXgivZ <- sqrt((shapeX / (2 * pi * x ^ 3))) * exp(- (shapeX * (x - meanX) ^ 2) / (2 * meanX ^ 2 * x))
    suppressWarnings(
      pXgivZ <- dinvGauss(x = x, lambda = shapeX, nu = meanX)
    )
    # -------------------------------------- Calculate
    # Check: shape and mean of inverse-Gaussian > 0 --
    pXgivZ[meanX <= 0] <- NA
    # -- Check: shape and mean of inverse-Gaussian > 0
  } else if (distX == "weibull") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- exp(eta_params[1])
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
    # pXgivZ <- (shapeX / scaleX) * ((x / scaleX) ^ (shapeX - 1)) * exp(- (x / scaleX) ^ shapeX)
    suppressWarnings(
      pXgivZ <- dweibull(x = x, shape = shapeX, scale = scaleX)
    )
    # -------------------------------------- Calculate
    # Check: shape and scale of Weibull > 0 ----------
    pXgivZ[scaleX <= 0] <- NA
    # ---------- Check: shape and scale of Weibull > 0
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
      # pXgivZ <- rateX * exp(- rateX * x)
      suppressWarnings(
        pXgivZ <- dexp(x = x, rate = rateX)
      )
    } else {
      # pXgivZ <- rateX ^ x * exp(- rateX) / factorial(x)
      suppressWarnings(
        pXgivZ <- dpois(x = x, lambda = rateX)
      )
    }
    # -------------------------------------- Calculate
    # Check: rate of exponential/Poisson > 0 ---------
    pXgivZ[rateX <= 0] <- NA
    # --------- Check: rate of exponential/Poisson > 0
  }
  return(pXgivZ)
}
