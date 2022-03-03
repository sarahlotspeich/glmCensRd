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
}

