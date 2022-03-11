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

part_deriv_pYgivXandZ <- function(y, x, z = NULL, distY, beta_params) {
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
    # Calculate P(Y|X,Z) -----------------------------
    eY <- as.numeric(y) - meanY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # ----------------------------- Calculate P(Y|X,Z)
    # Calculate partial derivatives ------------------
    all_d <- matrix(data = NA, nrow = nrow(pYgivXZ), ncol = length(beta_params))
    all_d[, 1] <- eY / (sigY ^ 2) * pYgivXZ # d/dbeta0
    all_d[, 2] <- x * all_d[, 1] # d/dbeta1
    # d/dbeta2, ..., d/dbetap
    if (!is.null(z)) {
      if (ncol(data.frame(z)) > 1) {
        for (c in 3:(length(beta_params) - 1)) {
          all_d[, c] <- data.frame(z)[, c] * all_d[, 1]
        }
      } else {
        all_d[, 3] <- z * all_d[, 1]
      }
    }
    all_d[, ncol(all_d)] <- (- (sigY ^ 2) + eY ^ 2) / (2 * (sigY ^ 2) ^ 2) * pYgivXZ # d/dsigma2
    # ------------------ Calculate partial derivatives
    # Give identifiable column names -----------------
    colnames(all_d) <- c(paste0("d_beta", 0:(length(beta_params) - 2)), "d_sigmaY2")
    # ----------------- Give identifiable column names
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
    eX <- data.matrix(x - meanX)
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
    # Give identifiable column names -----------------
    colnames(all_d) <- c(paste0("d_eta", 0:(length(eta_params) - 2)), "d_sigmaX2")
    # ----------------- Give identifiable column names
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
    beta_params <- params[1:(3 + length(Z))]
    # ------------------------------ Subset parameters
  } else if (distY == "binomial") {
    # Subset parameters ------------------------------
    beta_params <- params[1:(2 + length(Z))]
    # ------------------------------ Subset parameters
  }
  pYgivXZ <- calc_pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, beta_params = beta_params)
  d_pYgivXZ <- part_deriv_pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, beta_params = beta_params)
  d_loglik_theta <- d_pYgivXZ / matrix(data = pYgivXZ, nrow = length(pYgivXZ), ncol = length(beta_params))
  # Predictor model P(X|Z) ###########################
  # Subset parameters --------------------------------
  eta_params <- params[-c(1:length(beta_params))]
  # -------------------------------- Subset parameters
  pXgivZ <- calc_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)
  d_pXgivZ <- part_deriv_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)
  d_loglik_eta <- d_pXgivZ / matrix(data = pXgivZ, nrow = length(pXgivZ), ncol = length(eta_params))
  # Return deriv beta and eta side-by-side ----------
  d_loglik <- cbind(d_loglik_theta, d_loglik_eta)
  # ---------- Return deriv beta and eta side-by-side
  ####################################################
  # Derivatives of censored ##########################
  ####################################################
  # Integrate over joint P(Y,X,Z) --------------------
  dim_params <- length(c(beta_params, eta_params))
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
  integral_joint <- apply(X = cens_data, MARGIN = 1, FUN = integrate_joint_dens)
  integral_joint_wide <- matrix(data = rep(integral_joint, each = dim_params), ncol = dim_params, byrow = TRUE)
  # -------------------- Integrate over joint P(Y,X,Z)
  # Integrate over partial derivative ----------------
  part_deriv_joint_dens <- function(x, Yi, Zi, col = NULL) {
    ####################################################
    # Analysis model P(Y|X,Z) ##########################
    ####################################################
    pYgivXZ <- calc_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, beta_params = beta_params)
    d_pYgivXZ <- part_deriv_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, beta_params = beta_params)
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
    if (!is.null(col)) {
      return(d_pYXZ[, col])
    } else {
      return(d_pYXZ)
    }
  }
  integrate_part_deriv_joint_dens <- function(data_row) {
    data_row <- data.frame(t(data_row))
    return_mat <- matrix(data = NA, nrow = 1, ncol = dim_params)
    for (col in 1:ncol(return_mat)) {
      return_mat[, col] <- tryCatch(expr = integrate(f = part_deriv_joint_dens, lower = data_row[, W], upper = Inf, subdivisions = partX,
                                                     Yi = data_row[Y], Zi = data_row[, Z], col = col)$value,
                                    error = function(err) {0})
    }
    return(return_mat)
  }
  integral_d_loglik <- t(apply(X = cens_data, MARGIN = 1, FUN = integrate_part_deriv_joint_dens))
  # ---------------- Integrate over partial derivative
  # Divide integral of deriv by integral of P(Y,X,Z) -
  d_loglik <- rbind(d_loglik, integral_d_loglik / integral_joint_wide)
  # - Divide integral of deriv by integral of P(Y,X,Z)
  # Return matrix of derivatives ---------------------
  return(d_loglik)
  # --------------------- Return matrix of derivatives
}

# Calculate subject-specific derivatives -----------
## of the log-likelihood ---------------------------
## (Use as inputs into sandwich_B function) --------
d_theta <- part_deriv_loglik(params = params, Y = Y, W = W, D = D, Z = Z,
                             partX = partX, distY = distY, distX = distX, data = data)

sandwich_B <- function(d_theta) {
  # Construct the "meat" of the sandwich -------------
  B <- matrix(data = 0, nrow = ncol(d_theta), ncol = ncol(d_theta))
  for (c in 1:ncol(d_theta)) {
    for (r in c:ncol(d_theta)) {
      B[r, c] <- B[c, r] <- mean(d_theta[, c] * d_theta[, r])
    }
  }
  return(B)
}

second_deriv_pYgivXandZ <- function(y, x, z = NULL, distY, beta_params) {
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
    # Calculate P(Y|X,Z) -----------------------------
    eY <- as.numeric(y) - meanY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # ----------------------------- Calculate P(Y|X,Z)
    # Calculate partial derivatives ------------------
    d_pYgivXZ <- part_deriv_pYgivXandZ(y = y, x = x, z = z, distY = distY, beta_params = beta_params)
    # ------------------ Calculate partial derivatives
    # Calculate second derivatives -------------------
    num_deriv <- sum(length(beta_params) - 0:(length(beta_params) - 1))
    all_d2 <- matrix(data = NA, nrow = nrow(pYgivXZ), ncol = num_deriv)
    ## With respect to beta0 -------------------------
    start_col <- 1
    end_col <- 3 + ncol(data.frame(z))
    all_d2[, 1] <- (eY ^ 2 - sigY ^ 2) / ((sigY ^ 2) ^ 2) * pYgivXZ # d2/d2beta02
    all_d2[, 2] <- x * all_d2[, 1] # d2/dbeta0dbeta1
    # d/dbeta0dbeta2, ..., d/dbeta0dbetap
    if (!is.null(z)) {
      for (c in 1:ncol(data.frame(z))) {
        all_d2[, (2 + c)] <- data.frame(z)[, c] * all_d2[, 1]
      }
    }
    all_d2[, end_col] <- (- 2 - sigY ^ 2 + eY ^ 2) / (2 * sigY ^ 2) * d_pYgivXZ[, "d_beta0"] # d2/dbeta0dsigma2
    ## With respect to beta1 --------------------------
    start_col <- end_col + 1
    end_col <- end_col + (end_col - 1)
    all_d2[, start_col] <- x ^ 2 * all_d2[, 1] # d2/d2beta12
    # d/dbeta1dbeta2, ..., d/dbeta1dbetap
    if (!is.null(z)) {
      for (c in 1:ncol(data.frame(z))) {
        all_d2[, (start_col + c)] <- x * data.frame(z)[, c] * all_d2[, 1]
      }
    }
    all_d2[, end_col] <- (- 2 - sigY ^ 2 + eY ^ 2) / (2 * sigY ^ 2) * d_pYgivXZ[, "d_beta1"] # d2/dbeta1dsigma2
    ## With respect to beta21, ..., beta2p -----------
    if (!is.null(z)) {
      p <- ncol(data.frame(z))
      for (j in 1:p) {
        diff <- end_col - start_col
        start_col <- end_col + 1
        end_col <- end_col + diff
        all_d2[, start_col] <- data.frame(z)[, j] ^ 2 * all_d2[, 1] # d2/d2beta2j2
        if (p > 1) {
          for (c in 1:(p - j)) {
            all_d2[, (start_col + c)] <- data.frame(z)[, j] * data.frame(z)[, (c + 1)] * all_d2[, 1]
          }
        }
        all_d2[, end_col] <- (- 2 - sigY ^ 2 + eY ^ 2) / (2 * sigY ^ 2) * d_pYgivXZ[, (2 + j)] # d2/dbeta1dsigma2
      }
    }
    all_d2[, ncol(all_d2)] <- pYgivXZ / ((sigY ^ 2) ^ 2) * (((- sigY ^ 2 + eY ^ 2) / (2 * sigY ^ 2)) + 1) # d2/d2(sigma2)2
    # ------------------ Calculate partial derivatives
  }
  return(all_d2)
}

second_deriv_pXgivZ <- function(x, z = NULL, distX, eta_params) {
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
    # Calculate P(Y|X,Z) -----------------------------
    eX <- data.matrix(x - meanX)
    pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    # ----------------------------- Calculate P(Y|X,Z)
    # Calculate partial derivatives ------------------
    d_pXgivZ <- part_deriv_pXgivZ(x = x, z = z, distX = distX, eta_params = eta_params)
    # ------------------ Calculate partial derivatives
    # Calculate second derivatives -------------------
    num_deriv <- sum(length(eta_params) - 0:(length(eta_params) - 1))
    all_d2 <- matrix(data = NA, nrow = nrow(pXgivZ), ncol = num_deriv)
    ## With respect to eta0 --------------------------
    start_col <- 1
    end_col <- 2 + ncol(data.frame(z))
    all_d2[, 1] <- (eX ^ 2 - sigX ^ 2) / ((sigX ^ 2) ^ 2) * pXgivZ # d2/d2eta02
    # d/dbeta0deta1, ..., d/dbeta0dbetap
    if (!is.null(z)) {
      for (c in 1:ncol(data.frame(z))) {
        all_d2[, (1 + c)] <- data.frame(z)[, c] * all_d2[, 1]
      }
    }
    all_d2[, end_col] <- (- 2 - sigX ^ 2 + eX ^ 2) / (2 * sigX ^ 2) * d_pXgivZ[, 1] # d2/dbeta0dsigma2
    ## With respect to eta1, ..., etap ---------------
    if (!is.null(z)) {
      p <- ncol(data.frame(z))
      for (j in 1:p) {
        diff <- end_col - start_col
        start_col <- end_col + 1
        end_col <- end_col + diff
        all_d2[, start_col] <- data.frame(z)[, j] ^ 2 * all_d2[, 1] # d2/d2beta2j2
        if (p > 1) {
          for (c in 1:(p - j)) {
            all_d2[, (start_col + c)] <- data.frame(z)[, j] * data.frame(z)[, (c + 1)] * all_d2[, 1]
          }
        }
        all_d2[, end_col] <- (- 2 - sigX ^ 2 + eX ^ 2) / (2 * sigX ^ 2) * d_pXgivZ[, (1 + j)] # d2/detajdsigma2
      }
    }
    all_d2[, ncol(all_d2)] <- pXgivZ / ((sigX ^ 2) ^ 2) * (((- sigX ^ 2 + eX ^ 2) / (2 * sigX ^ 2)) + 1) # d2/d2(sigma2)2
    # ------------------ Calculate partial derivatives
  }
  return(all_d2)
}

second_deriv_loglik_wrt_beta <- function(pYgivXZ, d_pYgivXZ, d2_pYgivXZ, dimZ = 0) {
  # Create matrix to hold double derivatives of loglik -----------------
  d2_loglik_beta <- matrix(data = 0, nrow = nrow(d2_pYgivXZ), ncol = ncol(d2_pYgivXZ))
  ## With respect to beta0 /////////////////////////////////////////////
  start_col <- 1
  end_col <- 3 + dimZ
  # d2/d2beta02 of loglik ----------------------------------------------
  d2_loglik_beta[, 1] <- (d2_pYgivXZ[, start_col] * pYgivXZ - d_pYgivXZ[, 1] * d_pYgivXZ[, 1]) / (pYgivXZ ^ 2)
  # d2/d2beta12 of loglik ----------------------------------------------
  d2_loglik_beta[, (start_col + 1)] <- (d2_pYgivXZ[, (start_col + 1)] * pYgivXZ - d_pYgivXZ[, 1] * d_pYgivXZ[, 2]) / (pYgivXZ ^ 2)
  # d2/dbeta0dbeta2, ..., d2/dbeta0dbetap ------------------------------
  if (dimZ > 0) {
    for (c in 1:dimZ) {
      d2_loglik_beta[, (start_col + 1 + c)] <- (d2_pYgivXZ[, (start_col + 2 + c)] * pYgivXZ - d_pYgivXZ[, 1] * d_pYgivXZ[, (2 + c)]) / (pYgivXZ ^ 2)
    }
  }
  # d2/dbeta0dsigma2 of loglik -----------------------------------------
  d2_loglik_beta[, end_col] <- (d2_pYgivXZ[, end_col] * pYgivXZ - d_pYgivXZ[, 1] * d_pYgivXZ[, ncol(d_pYgivXZ)]) / (pYgivXZ ^ 2)
  ## ///////////////////////////////////////////// With respect to beta0
  ## With respect to beta1 /////////////////////////////////////////////
  diff <- end_col - start_col
  start_col <- end_col + 1
  end_col <- end_col + diff
  # d2/d2beta12 of loglik ----------------------------------------------
  d2_loglik_beta[, start_col] <- (d2_pYgivXZ[, start_col] * pYgivXZ - d_pYgivXZ[, 2] * d_pYgivXZ[, 2]) / (pYgivXZ ^ 2)
  # d2/dbeta1dbeta2, ..., d2/dbeta1dbetap ------------------------------
  if (dimZ > 0) {
    for (c in 1:dimZ) {
      d2_loglik_beta[, (start_col + c)] <- (d2_pYgivXZ[, (start_col + c)] * pYgivXZ - d_pYgivXZ[, 2] * d_pYgivXZ[, (2 + c)]) / (pYgivXZ ^ 2)
    }
  }
  # d2/dbeta0dsigma2 of loglik -----------------------------------------
  d2_loglik_beta[, end_col] <- (d2_pYgivXZ[, end_col] * pYgivXZ - d_pYgivXZ[, 2] * d_pYgivXZ[, ncol(d_pYgivXZ)]) / (pYgivXZ ^ 2)
  ## ///////////////////////////////////////////// With respect to beta1
  ## With respect to beta21, ..., beta2p ///////////////////////////////
  if (dimZ > 0) {
    for (j in 1:dimZ) {
      diff <- end_col - start_col
      start_col <- end_col + 1
      end_col <- end_col + diff
      # d2/d2betaj2 of loglik ----------------------------
      d2_loglik_beta[, start_col] <- (d2_pYgivXZ[, start_col] * pYgivXZ - d_pYgivXZ[, (2 + j)] * d_pYgivXZ[, (2 + j)]) / (pYgivXZ ^ 2)
      if (p > 1) {
        for (c in 1:(p - j)) {
          d2_loglik_beta[, (start_col + j)] <- (d2_pYgivXZ[, (start_col + j)] * pYgivXZ - d_pYgivXZ[, (2 + j)] * d_pYgivXZ[, (2 + j + c)]) / (pYgivXZ ^ 2)
        }
      }
      # d2/dbetajdsigma2 of loglik -----------------------
      d2_loglik_beta[, end_col] <- (d2_pYgivXZ[, end_col] * pYgivXZ - d_pYgivXZ[, (2 + j)] * d_pYgivXZ[, ncol(d_pYgivXZ)]) / (pYgivXZ ^ 2)
    }
  }
  ## /////////////////////////////// With respect to beta21, ..., beta2p
  ## With respect to sigma2 ////////////////////////////////////////////
  # d2/d2sigma22 of loglik ---------------------------
  d2_loglik_beta[, ncol(d2_loglik_beta)] <- (d2_pYgivXZ[, ncol(d2_pYgivXZ)] * pYgivXZ - d_pYgivXZ[, ncol(d_pYgivXZ)] * d_pYgivXZ[, ncol(d_pYgivXZ)]) / (pYgivXZ ^ 2)
  ## //////////////////////////////////////////// With respect to sigma2
  return(d2_loglik_beta)
}

second_deriv_loglik <- function(d_theta, params, Y, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data) {
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
    beta_params <- params[1:(3 + length(Z))]
    # ------------------------------ Subset parameters
  } else if (distY == "binomial") {
    # Subset parameters ------------------------------
    beta_params <- params[1:(2 + length(Z))]
    # ------------------------------ Subset parameters
  }
  pYgivXZ <- calc_pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, beta_params = beta_params)
  d_pYgivXZ <- d_theta[1:n1, 1:length(beta_params)]
  d2_pYgivXZ <- second_deriv_pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, beta_params = beta_params)
  d2_loglik_beta <- second_deriv_loglik_wrt_beta(pYgivXZ = pYgivXZ, d_pYgivXZ = d_pYgivXZ, d2_pYgivXZ = d2_pYgivXZ, dimZ = length(Z))
  # Predictor model P(X|Z) ###########################
  # Subset parameters --------------------------------
  eta_params <- params[-c(1:length(beta_params))]
  # -------------------------------- Subset parameters
  pXgivZ <- calc_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)
  d_pXgivZ <- d_theta[1:n1, - c(1:length(beta_params))]
  d2_pXgivZ <- second_deriv_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)
  ## Create matrix to hold double derivatives --------
  ## of loglik ---------------------------------------
  d2_loglik_eta <- matrix(data = 0, nrow = nrow(d2_pXgivZ), ncol = ncol(d2_pXgivZ))
  ## With respect to beta0 ---------------------------
  start_col <- 1
  end_col <- (2 + ncol(data.frame(z)))
  # d2/d2eta02 of loglik -----------------------------
  d2_loglik_eta[, start_col] <- (d2_pXgivZ[, start_col] * pXgivZ - d_pXgivZ[, 1] * d_pXgivZ[, 1]) / (pXgivZ ^ 2)
  if (!is.null(z)) {
    # d2/deta0deta1, ..., d2/deta0detap --------------
    for (c in 1:ncol(data.frame(z))) {
      d2_loglik_eta[, (start_col + c)] <- (d2_pXgivZ[, (start_col + c)] * pXgivZ - d_pXgivZ[, 1] * d_pXgivZ[, (1 + c)]) / (pXgivZ ^ 2)
    }
  }
  # d2/deta0dsigma2 of loglik ------------------------
  d2_loglik_eta[, end_col] <- (d2_pXgivZ[, end_col] * pXgivZ - d_pXgivZ[, 1] * d_pXgivZ[, ncol(d_pXgivZ)]) / (pXgivZ ^ 2)
  ## With respect to beta1 ---------------------------
  start_col <- (3 + c) + 1
  end_col <- 2 * (3 + c) - 1
  # d2/d2beta12 of loglik ----------------------------
  d2_loglik_beta[, start_col] <- (d2_pYgivXZ[, start_col] * pYgivXZ - d_pYgivXZ[, 2] * d_pYgivXZ[, 2]) / (pYgivXZ ^ 2)
  # d2/dbeta1dbeta2, ..., d2/dbeta1dbetap ------------
  if (!is.null(z)) {
    for (c in 1:ncol(data.frame(z))) {
      # d2/d2betaj2 of loglik
      d2_loglik_beta[, (start_col + c)] <- (d2_pYgivXZ[, (start_col + c)] * pYgivXZ - d_pYgivXZ[, 2] * d_pYgivXZ[, (2 + c)]) / (pYgivXZ ^ 2)
    }
    # d2/dbeta0dsigma2 of loglik
    d2_loglik_beta[, (start_col + c + 1)] <- (d2_pYgivXZ[, (2 + c)] * pYgivXZ - d_pYgivXZ[, 2] * d_pYgivXZ[, ncol(d_pYgivXZ)]) / (pYgivXZ ^ 2)
  }
  ## With respect to beta21, ..., beta2p -----------
  if (!is.null(z)) {
    p <- ncol(data.frame(z))
    for (j in 1:p) {
      diff <- end_col - start_col
      start_col <- end_col + 1
      end_col <- end_col + diff
      # d2/d2betaj2 of loglik ----------------------------
      d2_loglik_beta[, start_col] <- (d2_pYgivXZ[, start_col] * pYgivXZ - d_pYgivXZ[, (2 + j)] * d_pYgivXZ[, (2 + j)]) / (pYgivXZ ^ 2)
      if (p > 1) {
        for (c in 1:(p - j)) {
          d2_loglik_beta[, (start_col + j)] <- (d2_pYgivXZ[, (start_col + j)] * pYgivXZ - d_pYgivXZ[, (2 + j)] * d_pYgivXZ[, (2 + j + c)]) / (pYgivXZ ^ 2)
        }
      }
      # d2/dbetajdsigma2 of loglik -----------------------
      d2_loglik_beta[, end_col] <- (d2_pYgivXZ[, (start_col + c + 1)] * pYgivXZ - d_pYgivXZ[, (2 + j)] * d_pYgivXZ[, ncol(d_pYgivXZ)]) / (pYgivXZ ^ 2)
    }
  }
  # d2/d2sigma22 of loglik ---------------------------
  d2_loglik_beta[, ncol(d2_loglik_beta)] <- (d2_pYgivXZ[, ncol(d2_pYgivXZ)] * pYgivXZ - d_pYgivXZ[, ncol(d_pYgivXZ)] * d_pYgivXZ[, ncol(d_pYgivXZ)]) / (pYgivXZ ^ 2)

  # Return deriv beta and eta side-by-side ----------
  d_loglik <- cbind(d_loglik_beta, d_loglik_eta)
  # ---------- Return deriv beta and eta side-by-side
  ####################################################
  # Derivatives of censored ##########################
  ####################################################
  # Integrate over joint P(Y,X,Z) --------------------
  dim_params <- length(c(beta_params, eta_params))
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
  integral_joint <- apply(X = cens_data, MARGIN = 1, FUN = integrate_joint_dens)
  integral_joint_wide <- matrix(data = rep(integral_joint, each = dim_params), ncol = dim_params, byrow = TRUE)
  # -------------------- Integrate over joint P(Y,X,Z)
  # Integrate over partial derivative ----------------
  part_deriv_joint_dens <- function(x, Yi, Zi, col = NULL) {
    ####################################################
    # Analysis model P(Y|X,Z) ##########################
    ####################################################
    pYgivXZ <- calc_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, beta_params = beta_params)
    d_pYgivXZ <- part_deriv_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, beta_params = beta_params)
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
    if (!is.null(col)) {
      return(d_pYXZ[, col])
    } else {
      return(d_pYXZ)
    }
  }
  integrate_part_deriv_joint_dens <- function(data_row) {
    data_row <- data.frame(t(data_row))
    return_mat <- matrix(data = NA, nrow = 1, ncol = dim_params)
    for (col in 1:ncol(return_mat)) {
      return_mat[, col] <- tryCatch(expr = integrate(f = part_deriv_joint_dens, lower = data_row[, W], upper = Inf, subdivisions = partX,
                                                     Yi = data_row[Y], Zi = data_row[, Z], col = col)$value,
                                    error = function(err) {0})
    }
    return(return_mat)
  }
  integral_d_loglik <- t(apply(X = cens_data, MARGIN = 1, FUN = integrate_part_deriv_joint_dens))
  # ---------------- Integrate over partial derivative
  # Divide integral of deriv by integral of P(Y,X,Z) -
  d_loglik <- rbind(d_loglik, integral_d_loglik / integral_joint_wide)
  # - Divide integral of deriv by integral of P(Y,X,Z)
  # Return matrix of derivatives ---------------------
  return(d_loglik)
  # --------------------- Return matrix of derivatives
}
