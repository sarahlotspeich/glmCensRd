calc_pYgivXandZ <- function(y, x, z = NULL, distY, theta_params) {
  if (distY == "normal") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- theta_params[1] + theta_params[2] * matrix(data = x, ncol = 1)
    if (!is.null(z)) {
      beta2 <- theta_params[-c(1:2, length(theta_params))]
      if (length(beta2) == 1) {
        meanY <- meanY + beta2 * z
      } else {
        meanY <- meanY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
      }
    }
    ## Estimate sqrt(variance) directly --------------
    sigY <- theta_params[length(theta_params)]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    eY <- as.numeric(y) - meanY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # -------------------------------------- Calculate
  } else if (distY == "binomial") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- theta_params[1] + theta_params[2] * matrix(data = x, ncol = 1)
    if (!is.null(z)) {
      beta2 <- theta_params[-c(1:2)]
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

part_deriv_pYgivXandZ <- function(y, x, z = NULL, distY, theta_params) {
  if (distY == "normal") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- theta_params[1] + theta_params[2] * matrix(data = x, ncol = 1)
    if (!is.null(z)) {
      beta2 <- theta_params[-c(1:2, length(theta_params))]
      if (length(beta2) == 1) {
        meanY <- meanY + beta2 * z
      } else {
        meanY <- meanY + as.numeric(data.matrix(z) %*% matrix(data = beta2, ncol = 1))
      }
    }
    ## Estimate sqrt(variance) directly --------------
    sigY <- theta_params[length(theta_params)]
    # --------------------------------- Get parameters
    # Calculate P(Y|X,Z) -----------------------------
    eY <- as.numeric(y) - meanY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # ----------------------------- Calculate P(Y|X,Z)
    # Calculate partial derivatives ------------------
    all_d <- matrix(data = NA, nrow = nrow(pYgivXZ), ncol = length(theta_params))
    all_d[, 1] <- eY / (sigY ^ 2) * pYgivXZ # d/dbeta0
    all_d[, 2] <- x * all_d[, 1] # d/dbeta1
    # d/dbeta2, ..., d/dbetap
    if (!is.null(z)) {
      if (ncol(data.frame(z)) > 1) {
        for (c in 3:(length(theta_params) - 1)) {
          all_d[, c] <- data.frame(z)[, c] * all_d[, 1]
        }
      } else {
        all_d[, 3] <- z * all_d[, 1]
      }
    }
    all_d[, ncol(all_d)] <- (- (sigY ^ 2) + eY ^ 2) / (2 * (sigY ^ 2) ^ 2) * pYgivXZ # d/dsigma2
    # ------------------ Calculate partial derivatives
  } else if (distY == "binomial") {
    # Get parameters ---------------------------------
    ## Construct mean --------------------------------
    meanY <- theta_params[1] + theta_params[2] * matrix(data = x, ncol = 1)
    if (!is.null(z)) {
      beta2 <- theta_params[-c(1:2)]
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
  return(all_d)
}

part_deriv_pXgivZ <- function(x, z = NULL, distX, eta_params) {
  if (distX == "normal") {
    # Calculate P(Y|X,Z) -----------------------------
    pXgivZ <- calc_pXgivZ(x = x, z = z, distX = distX, eta_params = eta_params)
    # ----------------------------- Calculate P(Y|X,Z)
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
    eX <- x - meanX
    pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    # --------------------------------- Get parameters
    # Calculate partial derivatives ------------------
    all_d <- matrix(data = NA, nrow = nrow(pXgivZ), ncol = length(eta_params))
    all_d[, 1] <- eX / (sigX ^ 2) * pXgivZ # d/deta0
    #all_d[, 2] <- x * all_d[, 1] # d/deta1
    # d/deta1, ..., d/detap
    if (!is.null(z)) {
      if (ncol(data.frame(z)) > 1) {
        for (c in 2:(length(eta_params) - 1)) {
          all_d[, c] <- data.frame(z)[, c] * all_d[, 1]
        }
      } else {
        all_d[, 2] <- z * all_d[, 1]
      }
    }
    all_d[, ncol(all_d)] <- (- (sigX ^ 2) + eX ^ 2) / (2 * (sigX ^ 2) ^ 2) * pXgivZ # d/dsigma2
    # ------------------ Calculate partial derivatives
  }
  return(all_d)
}

part_deriv_loglik <- function(params, Y, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data) {
  ####################################################
  # Pre-processing ###################################
  ####################################################
  # Subset data to relevant, user-specified columns --
  data <- data[, c(Y, W, D, Z)]
  # Create variable X = W ----------------------------
  data <- cbind(data, X = data[, W])
  ## Make it NA for censored (D = 0) -----------------
  data[data[, D] == 0, "X"] <- NA
  ## < define predictor column name > ----------------
  X <- "X"
  ## ---------------- < define predictor column name >
  # < number of uncensored subjects > ----------------
  n1 <- sum(data[, D]) # -----------------------------
  # ---------------- < number of uncensored subjects >
  # Reordered data to be uncensored first ------------
  data <- data[order(data[, D], decreasing = TRUE), ]
  # ------------ Reordered data to be uncensored first
  # Create subset of uncensored subjects' data -------
  uncens_data <- data[1:n1, ]
  # ------- Create subset of uncensored subjects' data
  # Create subset of censored subjects' data ---------
  cens_data <- data[-c(1:n1), ]
  # --------- Create subset of censored subjects' data
  ####################################################
  # Derivatives of uncensored ########################
  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  if (distY == "normal") {
    # Subset parameters ------------------------------
    theta_params <- params[1:(3 + length(Z))]
    # ------------------------------ Subset parameters
  } else if (distY == "binomial") {
    # Subset parameters ------------------------------
    theta_params <- params[1:(2 + length(Z))]
    # ------------------------------ Subset parameters
  }
  pYgivXZ <- calc_pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, theta_params = theta_params)
  d_pYgivXZ <- part_deriv_pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, theta_params = theta_params)
  d_loglik_theta <- d_pYgivXZ / matrix(data = pYgivXZ, nrow = length(pYgivXZ), ncol = length(theta_params))
  # Predictor model P(X|Z) ###########################
  # Subset parameters --------------------------------
  eta_params <- params[-c(1:length(theta_params))]
  # -------------------------------- Subset parameters
  pXgivZ <- calc_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)
  d_pXgivZ <- part_deriv_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)
  d_loglik_eta <- d_pXgivZ / matrix(data = pXgivZ, nrow = length(pXgivZ), ncol = length(eta_params))
  ####################################################
  # Derivatives of censored ##########################
  ####################################################
  # Integrate over joint P(Y,X,Z) --------------------
  joint_dens <- function(x, Yi, Zi) {
    ####################################################
    # Analysis model P(Y|X,Z) ##########################
    ####################################################
    pYgivXZ <- calc_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, theta_params = theta_params)

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
  # -------------------- Integrate over joint P(Y,X,Z)
  # Integrate over partial derivative ----------------
  part_deriv_joint_dens <- function(x, Yi, Zi, col) {
    ####################################################
    # Analysis model P(Y|X,Z) ##########################
    ####################################################
    pYgivXZ <- calc_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, theta_params = theta_params)
    d_pYgivXZ <- part_deriv_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, theta_params = theta_params)
    ####################################################
    # Predictor model P(X|Z) ###########################
    ####################################################
    pXgivZ <- calc_pXgivZ(x = x, z = Zi, distX = distX, eta_params = eta_params)
    d_pXgivZ <- part_deriv_pXgivZ(x = x, z = Zi, distX = distX, eta_params = eta_params)
    ####################################################
    # Partial derivatives ##############################
    ####################################################
    d_pYXZ <- cbind(d_pYgivXZ * matrix(data = pXgivZ, nrow = length(pXgivZ), ncol = ncol(d_pYgivXZ), byrow = FALSE),
                    matrix(data = pYgivXZ, nrow = length(pYgivXZ), ncol = ncol(d_pXgivZ), byrow = FALSE) * d_pXgivZ)
    return(d_pYXZ[, col])
  }
  integrate_part_deriv_joint_dens <- function(data_row) {
    data_row <- data.frame(t(data_row))
    return_mat <- matrix(data = NA, nrow = 1, ncol = c(theta_params, eta_params))
    for (col in 1:ncol(return_mat)) {
      return_mat[, col] <- tryCatch(expr = integrate(f = part_deriv_joint_dens, lower = data_row[, W], upper = Inf, subdivisions = partX,
                                                     Yi = data_row[Y], Zi = data_row[, Z], col = col)$value,
                                    error = function(err) {0})
    }
    return(return_mat)
  }

  integrate_part_deriv_joint_dens(cens_data[1, ])

  integrate_d_loglik <- matrix(data = NA, nrow = nrow(cens_data), ncol = length(c(theta_params, eta_params)))
  # ---------------- Integrate over partial derivative
  # Return matrix of derivatives ---------------------
  return(cbind(d_loglik_theta, d_loglik_eta))
  # --------------------- Return matrix of derivatives
}

sandwich_B <- function(y, x, z = NULL, distY, theta_params, distX, eta_params) {
  d_theta <- part_deriv_pYgivXandZ(y = y, x = x, z = z, distY = distY, theta_params = theta_params)
  d_eta <- part_deriv_pXgivZ(x = x, z = z, distX = distX, eta_params = eta_params)
  d_theta_eta <- cbind(d_theta, d_eta)
  B <- matrix(data = 0, nrow = ncol(d_theta_eta), ncol = ncol(d_theta_eta))
  for (c in 1:ncol(d_theta_eta)) {
    for (r in c:ncol(d_theta_eta)) {
      B[r, c] <- B[c, r] <- mean(d_theta_eta[, c] * d_theta_eta[, r])
    }
  }
  return(B)
}
