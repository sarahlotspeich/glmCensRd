pYgivXandZ <- function(y, x, z = NULL, distY, theta_params) {
  if (distY == "normal") {
    # Get parameters ---------------------------------
    beta0 <- theta_params[1]
    beta1 <- theta_params[2]
    if (!is.null(z)) { beta2 <- theta_params[3:(2 + length(z))] }
    sigY <- theta_params[(2 + length(Z)) + 1]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muY <- beta0 + beta1 * matrix(data = x, ncol = 1)
    if (length(z) > 0) {
      muY <- muY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
    }
    muY <- data.matrix(muY)
    eY <- as.numeric(y) - muY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # -------------------------------------- Calculate
  } else if (distY == "binomial") {
    # Get parameters ---------------------------------
    beta0 <- theta_params[1]
    beta1 <- theta_params[2]
    if (!is.null(z)) { beta2 <- theta_params[3:(2 + length(z))] }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muY <- beta0 + beta1 * matrix(data = x, ncol = 1)
    if (length(z) > 0) {
      muY <- muY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
    }
    muY <- data.matrix(muY)
    pYgivXZ <- exp(- (1 - y) * muY) / (1 + exp(- muY))
    # -------------------------------------- Calculate
  }
  return(pYgivXZ)
}

pXgivZ <- function(x, z = NULL, distX, eta_params) {
  if (distX == "normal") {
    # Get parameters ---------------------------------
    eta0 <- eta_params[1]
    if (!is.null(z)) { eta1 <- eta_params[2:(1 + length(z))] }
    sigX <- eta_params[(1 + length(Z)) + 1]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muX <- eta0
    if (length(z) > 0) {
      muX <- muX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    eX <- x - muX
    pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    # -------------------------------------- Calculate
  } else if (distX == "log-normal") {
    # Get parameters ---------------------------------
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(z))] }
    sigX <- eta_params[(1 + length(Z)) + 1]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muX <- eta0
    if (length(z) > 0) {
      muX <- muX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    pXgivZ <- (1 / (x * sigX * sqrt(2 * pi))) * exp(- (log(x) - muX) ^ 2 / (2 * sigX ^ 2))
    # -------------------------------------- Calculate
  } else if (distX == "gamma") {
    # Get parameters ---------------------------------
    ## Estimate shape directly -----------------------
    shapeX <- eta_params[1]
    #eta0 <- eta_params[1]
    #if (!is.null(z)) { eta1 <- eta_params[2:(1 + length(z))] }
    #shapeX <- eta0
    #if (length(z) > 0) {
    #  shapeX <- shapeX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    #}
    #eta_params <- eta_params[-c(1:(1 + length(z)))]
    ## Construct mean --------------------------------
    eta0 <- eta_params[2]
    if (!is.null(z)) { eta1 <- eta_params[3:(2 + length(z))] }
    muX <- eta0
    if (length(z) > 0) {
      muX <- muX + as.numeric(data.matrix(z) %*% matrix(data = eta1, ncol = 1))
    }
    ## Construct scale -------------------------------
    scaleX <- muX / shapeX
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- (1 / gamma(shapeX)) * scaleX ^ (- shapeX) * (x ^ (shapeX - 1)) * exp(- x / scaleX)
    # -------------------------------------- Calculate
  } else if (distX == "inverse-gaussian") {
    # Get parameters ---------------------------------
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(Z))] }
    shapeX <- exp(eta_params[(1 + length(Z)) + 1])
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muX <- eta0
    if (length(Z) > 0) {
      muX <- exp(muX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1)))
    }
    pXgivZ <- sqrt((shapeX / (2 * pi * uncens_data[, X]^3))) * exp(-1 * (shapeX * (uncens_data[, X] - muX)^2) / (2 * muX^2 * uncens_data[, X]))
    # -------------------------------------- Calculate
  } else if (distX == "weibull") {
    # Get parameters ---------------------------------
    shapex <- eta_params[1]
    eta_params <- eta_params[-1]
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(Z))] }
    scalex <- eta0
    if (length(Z) > 0) {
      scalex <- scalex + data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1)
    }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    pXgivZ <- (shapex / scalex) * ((uncens_data[, X] / scalex) ^ (shapex - 1)) * exp(-1 * (uncens_data[, X] / scalex)^ shapex)
    # -------------------------------------- Calculate
  } else if (distX == "exponential") {
    # Get parameters ---------------------------------
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(Z))] }
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    lambdaX <- eta0
    if (length(Z) > 0) {
      lambdaX <- lambdaX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
    }
    pXgivZ <- lambdaX * exp(- lambdaX * uncens_data[, X])
    # -------------------------------------- Calculate
  } else if (distX == "poisson") {
    # Get parameters ---------------------------------
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(Z))] }
    # --------------------------------- Get parameters
    muX <- eta0
    if (length(Z) > 0) {
      muX <- muX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
    }
    pXgivZ <- ((muX^uncens_data[, X]) * (exp(-1 * muX))) / factorial(uncens_data[, X])
    # -------------------------------------- Calculate
  }
  return(pXgivZ)
}