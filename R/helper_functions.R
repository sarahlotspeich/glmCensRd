pYgivXandZ <- function(y, x, z = NULL, distY, theta_params) {
  if (distY == "normal") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- theta_params[1] + theta_params[2] * matrix(data = x, ncol = 1)
    if (length(z) > 0) {
      beta2 <- theta_params[3:(2 + length(z))]
      meanY <- meanY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
    }
    meanY <- data.matrix(muY)
    ## Estimate sqrt(variance) directly --------------
    sigY <- theta_params[(2 + length(Z)) + 1]
    # --------------------------------- Get parameters
    eY <- as.numeric(y) - meanY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # -------------------------------------- Calculate
  } else if (distY == "binomial") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- theta_params[1] + theta_params[2] * matrix(data = x, ncol = 1)
    if (length(z) > 0) {
      beta2 <- theta_params[3:(2 + length(z))]
      meanY <- meanY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
    }
    meanY <- data.matrix(meanY)
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pYgivXZ <- exp(- (1 - y) * meanY) / (1 + exp(meanY))
    # -------------------------------------- Calculate
  }
  return(pYgivXZ)
}

pXgivZ <- function(x, z = NULL, distX, eta_params) {
  if (distX == "normal") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanX <- eta_params[1]
    if (length(z) > 0) {
      eta1 <- eta_params[2:(1 + length(z))]
      meanX <- meanX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    ## Estimate sqrt(variance) directly --------------
    sigX <- eta_params[(1 + length(Z)) + 1]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    eX <- x - meanX
    pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    # -------------------------------------- Calculate
  } else if (distX == "log-normal") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanX <- eta_params[1]
    if (length(z) > 0) {
      eta1 <- eta_params[2:(1 + length(z))]
      meanX <- meanX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    ## Estimate sqrt(variance) directly --------------
    sigX <- eta_params[(1 + length(Z)) + 1]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- (1 / (x * sigX * sqrt(2 * pi))) * exp(- (log(x) - meanX) ^ 2 / (2 * sigX ^ 2))
    # -------------------------------------- Calculate
  } else if (distX == "gamma") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- eta_params[1]
    ## Construct mean --------------------------------
    meanX <- eta_params[2]
    if (length(z) > 0) {
      eta1 <- eta_params[3:(2 + length(z))]
      meanX <- meanX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
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
    if (length(z) > 0) {
      meanX <- meanX + as.numeric(data.matrix(z) %*% matrix(data = eta_params[3:(2 + length(z))], ncol = 1))
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
    if (length(Z) > 0) {
      eta1 <- eta_params[3:(2 + length(z))]
      scaleX <- scaleX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- (shapeX / scaleX) * ((x / scaleX) ^ (shapeX - 1)) * exp(- (x / scaleX) ^ shapeX)
    # -------------------------------------- Calculate
  } else if (distX == "exponential") {
    # Get parameters ---------------------------------
    ## Construct rate  -------------------------------
    rateX <-eta_params[1]
    if (length(z) > 0) {
      eta1 <- eta_params[2:(1 + length(z))]
      rateX <- rateX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- rateX * exp(- rateX * x)
    # -------------------------------------- Calculate
  } else if (distX == "poisson") {
    # Get parameters ---------------------------------
    ## Construct rate  -------------------------------
    rateX <-eta_params[1]
    if (length(z) > 0) {
      eta1 <- eta_params[2:(1 + length(z))]
      rateX <- rateX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- ((rateX ^ x) * (exp(- rateX))) / factorial(x)
    # -------------------------------------- Calculate
  }
  return(pXgivZ)
}