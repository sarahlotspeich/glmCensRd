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
    pYgivXZ <- exp(- (1 - y) * meanY) / (1 + exp(meanY))
    # -------------------------------------- Calculate
  }
  return(pYgivXZ)
}

calc_pXgivZ <- function(x, z = NULL, distX, eta_params) {
  if (distX == "normal") {
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
    pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    # -------------------------------------- Calculate
  } else if (distX == "log-normal") {
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
    pXgivZ <- (1 / (x * sigX * sqrt(2 * pi))) * exp(- (log(x) - meanX) ^ 2 / (2 * sigX ^ 2))
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
  } else if (distX == "exponential") {
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
    pXgivZ <- rateX * exp(- rateX * x)
    # -------------------------------------- Calculate
  } else if (distX == "poisson") {
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
    pXgivZ <- rateX ^ x * exp(- rateX) / factorial(x)
    # -------------------------------------- Calculate
  }
  return(pXgivZ)
}

calc_indiv_loglik <- function(params, Y, X, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data) {
  ####################################################
  # Pre-processing ###################################
  ####################################################
  # < number of uncensored subjects > ----------------
  n1 <- sum(data[, D]) # -----------------------------
  # ---------------- < number of uncensored subjects >
  # Reordered data to be uncensored first ------------
  data <- data[order(data[, D], decreasing = TRUE), ]
  # ------------ Reordered data to be uncensored first
  # Create subset of uncensored subjects' data -------
  uncens_data <- data[1:n1, ]
  # ------- Create subset of uncensored subjects' data
  # Create subset of censored subjects' data -------
  cens_data <- data[-c(1:n1), ]
  # ------- Create subset of censored subjects' data
  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  ####################################################
  if (distY == "normal") {
    # Subset parameters ------------------------------
    beta_params <- params[1:(3 + length(Z))]
    # ------------------------------ Subset parameters
  } else if (distY == "binomial") {
    # Subset parameters ------------------------------
    beta_params <- params[1:(2 + length(Z))]
    # ------------------------------ Subset parameters
  }
  pYgivXZ <- calc_pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, beta_params = beta_params)

  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  # Subset parameters --------------------------------
  eta_params <- params[-c(1:length(beta_params))]
  # -------------------------------- Subset parameters
  # Check for parameters outside domain --------------
  if (distX == "gamma") {
    # Shape and scale of Gamma > 0 -------------------
    shapeX <- eta_params[1]
    meanX <- eta_params[2]
    if (length(Z) > 0) {
      meanX <- meanX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta_params[3:(2 + length(Z))], ncol = 1))
    }
    scaleX <- meanX / shapeX
    if (any(c(shapeX, scaleX) <= 0)) { return(99999)}
    # ------------------- Shape and scale of Gamma > 0
  } else if (distX == "inverse-gaussian") {
    # Shape and mean of inverse-Gaussian both > 0 ----
    shapeX <- eta_params[1]
    meanX <- eta_params[2]
    if (length(Z) > 0) {
      meanX <- meanX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta_params[3:(2 + length(Z))], ncol = 1))
    }
    if (any(c(shapeX, meanX) <= 0)) { return(99999)}
    # --------------------------------- Get parameters
  } else if (distX == "weibull") {
    # Shape and scale of Weibull both > 0 ------------
    shapeX <- eta_params[1]
    scaleX <- eta_params[2]
    if (length(Z) > 0) {
      eta1 <- eta_params[3:(2 + length(Z))]
      scaleX <- scaleX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
    }
    if (any(c(shapeX, scaleX) <= 0)) { return(99999)}
  } else if (distX %in% c("exponential", "poisson")) {
    # Rate of Exponential or Poisson > 0 -------------
    rateX <- eta_params[1]
    if (length(Z) > 0) {
      eta1 <- eta_params[2:(1 + length(Z))]
      if (length(eta1) == 1) {
        rateX <- rateX + eta1 * uncens_data[, Z]
      } else {
        rateX <- rateX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
      }
    }
    if (any(rateX <= 0)) { return (99999) }
  }
  pXgivZ <- calc_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)

  ####################################################
  # Calculate joint density P(Y,X,Z) #################
  ####################################################
  uncens_data <- cbind(uncens_data, jointP = 1)
  uncens_data[, "jointP"] <- pYgivXZ * pXgivZ
  #uncens_data <- data.frame(cbind(uncens_data, jointP = pYgivXZ * pXgivZ))

  ####################################################
  # Calculate the log-likelihood #####################
  ####################################################
  # Log-likelihood contribution of uncensored X ------
  ll <- log(uncens_data[, "jointP"])
  # ------ Log-likelihood contribution of uncensored X
  # Log-likelihood contribution of censored X --------
  joint_dens <- function(x, Yi, Zi) {
    ####################################################
    # Analysis model P(Y|X,Z) ##########################
    ####################################################
    pYgivXZ <- calc_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, beta_params = beta_params)

    ####################################################
    # Predictor model P(X|Z) ###########################
    ####################################################
    pXgivZ <- calc_pXgivZ(x = x, z = Zi, distX = distX, eta_params = eta_params)

    ####################################################
    # Joint density P(Y,X,Z) ###########################
    ####################################################
    return(pYgivXZ * pXgivZ)
  }
  integrate_joint_dens <- function(data_row) {
    data_row <- data.frame(t(data_row))
    return(
      tryCatch(expr = integrate(f = joint_dens, lower = data_row[, W], upper = Inf, subdivisions = partX,
                                Yi = data_row[Y], Zi = data_row[, Z])$value,
               error = function(err) {0})
    )
  }
  integral <- apply(X = cens_data, MARGIN = 1, FUN = integrate_joint_dens)
  log_integral <- log(integral)
  log_integral[log_integral == -Inf] <- 0
  ll <- append(ll, log_integral)
  # -------- Log-likelihood contribution of censored X
  # Return (-1) x log-likelihood for use with nlm() --
  ll[ll == - Inf] <- 0
  return(ll)
  # -- Return (-1) x log-likelihood for use with nlm()
}

calc_deriv_loglik <- function(mle, Y, X, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data, j = NULL) {
  p <- length(mle)
  eps <- mle * (10 ^ (- 4))

  if (is.null(j)) {
    # Create matrix to save derivatives in
    d_theta <- matrix(data = 0, nrow = nrow(data), ncol = p, byrow = FALSE)
    for (j in 1:p) {
      # Create jth euclidean vector
      ej <- matrix(data = 0, nrow = p, ncol = 1)
      ej[j] <- 1

      # MLE - eps
      mle0 <- matrix(data = mle - eps * ej, nrow = p, ncol = 1)
      l0 <- calc_indiv_loglik(params = mle0, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)

      # MLE + eps
      mle1 <- matrix(data = mle + eps * ej, nrow = p, ncol = 1)
      l1 <- calc_indiv_loglik(params = mle1, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)

      # Estimate deriv = (l(beta + eps) - l(beta - eps)) / (2 * eps)
      d_theta[, j] <- (l1 - l0) / (2 * eps[j])
    }
  } else {
    # Create matrix to save derivatives in
    d_theta <- matrix(data = 0, nrow = nrow(data), ncol = 1, byrow = FALSE)

    # Create jth euclidean vector
    ej <- matrix(data = 0, nrow = p, ncol = 1)
    ej[j] <- 1

    # MLE - eps
    mle0 <- matrix(data = mle - eps * ej, nrow = p, ncol = 1)
    l0 <- calc_indiv_loglik(params = mle0, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)

    # MLE + eps
    mle1 <- matrix(data = mle + eps * ej, nrow = p, ncol = 1)
    l1 <- calc_indiv_loglik(params = mle1, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)

    # Estimate deriv = (l(beta + eps) - l(beta - eps)) / (2 * eps)
    d_theta[, j] <- (l1 - l0) / (2 * eps[j])
  }
  return(d_theta)
}

calc_deriv_loglik <- function(mle, Y, X, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data, j = NULL) {
  p <- length(mle)
  hn <- nrow(data) ^ (- 1 / 2)
  # Calculate the log-likelihood contributions at the MLE
  l <- calc_indiv_loglik(params = mle, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
  if (is.null(j)) {
    # Create matrix to save
    d_theta <- matrix(data = - l, nrow = nrow(data), ncol = p, byrow = FALSE)
    for (j in 1:p) {
      ej <- matrix(data = 0, nrow = p, ncol = 1)
      ej[j] <- 1
      mle_ <- matrix(data = mle + hn * ej, nrow = p, ncol = 1)
      l_ <- calc_indiv_loglik(params = mle_, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
      d_theta[, j] <- d_theta[, j] + l_
    }
    d_theta <- (1 / hn) * d_theta
  } else {
    # Create matrix to save
    d_theta <- matrix(data = - l, nrow = nrow(data), ncol = 1, byrow = FALSE)
    ej <- matrix(data = 0, nrow = p, ncol = 1)
    ej[j] <- 1
    mle_ <- matrix(data = mle + hn * ej, nrow = p, ncol = 1)
    l_ <- calc_indiv_loglik(params = mle_, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
    d_theta[, j] <- d_theta[, j] + l_
    d_theta <- (1 / hn) * d_theta
  }
  return(d_theta)
}

calc_deriv2_loglik <- function(mle, Y, X, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data, j = NULL, k = NULL) {
  p <- length(mle)
  hn <- nrow(data) ^ (- 1 / 2)
  # Calculate the log-likelihood contributions at the MLE
  l <- calc_indiv_loglik(params = mle, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
  if (is.null(j)) {
    if (!is.null(k)) {
      # Create matrix to save
      d_theta <- matrix(data = l, nrow = nrow(data), ncol = p, byrow = FALSE)
      # Perturb kth element of mle
      ek <- matrix(data = 0, nrow = p, ncol = 1)
      ek[k] <- 1
      mle_k <- matrix(data = mle + hn * ek, nrow = p, ncol = 1)
      l_k <- calc_indiv_loglik(params = mle_k, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
      for (j in 1:p) {
        ej <- matrix(data = 0, nrow = p, ncol = 1)
        ej[j] <- 1
        mle_j <- matrix(data = mle + hn * ej, nrow = p, ncol = 1)
        mle_jk <- matrix(data = mle + hn * ej + hn * ek, nrow = p, ncol = 1)
        l_j <- calc_indiv_loglik(params = mle_j, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
        l_jk <- calc_indiv_loglik(params = mle_jk, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
        d_theta[, j] <- d_theta[, j] + l_jk - l_j - l_k
      }
      d_theta <- (1 / (hn ^ 2)) * d_theta
    } else {
      for (k in 1:p) {
        # Create matrix to save
        d_theta_k <- matrix(data = l, nrow = nrow(data), ncol = p, byrow = FALSE)
        # Perturb kth element of mle
        ek <- matrix(data = 0, nrow = p, ncol = 1)
        ek[k] <- 1
        mle_k <- matrix(data = mle + hn * ek, nrow = p, ncol = 1)
        l_k <- calc_indiv_loglik(params = mle_k, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
        for (j in 1:p) {
          ej <- matrix(data = 0, nrow = p, ncol = 1)
          ej[j] <- 1
          mle_j <- matrix(data = mle + hn * ej, nrow = p, ncol = 1)
          mle_jk <- matrix(data = mle + hn * ej + hn * ek, nrow = p, ncol = 1)
          l_j <- calc_indiv_loglik(params = mle_j, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
          l_jk <- calc_indiv_loglik(params = mle_jk, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
          d_theta_k[, j] <- d_theta_k[, j] + l_jk - l_j - l_k
        }
        if (k > 1) {
          d_theta <- cbind(d_theta, d_theta_k)
        } else {
          d_theta <- d_theta_k
        }
      }
      d_theta <- (1 / (hn ^ 2)) * d_theta
    }
  } else {
    # Create matrix to save
    d_theta <- matrix(data = l, nrow = nrow(data), ncol = 1, byrow = FALSE)
    ej <- ek <- matrix(data = 0, nrow = p, ncol = 1)
    ej[j] <- 1
    ek[k] <- 1
    mle_j <- matrix(data = mle + hn * ej, nrow = p, ncol = 1)
    mle_k <- matrix(data = mle + hn * ek, nrow = p, ncol = 1)
    mle_jk <- matrix(data = mle + hn * ej + hn * ek, nrow = p, ncol = 1)
    l_j <- calc_indiv_loglik(params = mle_j, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
    l_k <- calc_indiv_loglik(params = mle_k, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
    l_jk <- calc_indiv_loglik(params = mle_jk, Y = Y, X = X, W = W, D = D, Z = Z, partX = partX, distY = distX, distX = distX, data = data)
    d_theta[, j] <- d_theta[, j] + l_jk - l_j - l_k
    d_theta <- (1 / (hn ^ 2)) * d_theta
  }
  return(d_theta)
}
