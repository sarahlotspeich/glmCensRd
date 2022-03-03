#' Observed-data log-likelihood

#' @param params Parameter values.
#' @param Y Name of outcome variable.
#' @param X Name of censored predictor variable.
#' @param W Name of observed (i.e., censored) version of \code{X}.
#' @param D Name of event indicator, defined to be = 1 if \code{X} was uncensored.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{NULL}
#' @param partX Size of partition of unobserved \code{X} for censored subjects. Default is \code{50}.
#' @param distY Distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}.
#' @param distX Distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}.
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{W}, \code{D}, and \code{Z}.

#' @return A vector containing the approximate integral over the joint density of censored subjects from \code{complete_data_cens}.
#'
loglik <- function(params, Y, X, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data) {
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
    theta_params <- params[1:(3 + length(Z))]
    # ------------------------------ Subset parameters
  } else if (distY == "binomial") {
    # Subset parameters ------------------------------
    theta_params <- params[1:(2 + length(Z))]
    # ------------------------------ Subset parameters
  }
  pYgivXZ <- pYgivXandZ(y = uncens_data[, Y], x = uncens_data[, X], z = uncens_data[, Z], distY = distY, theta_params = theta_params)
  
  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  if (distX == "normal") {
    # Get parameters ---------------------------------
    eta_params <- params[-c(1:length(theta_params))]
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(Z))] }
    sigX <- eta_params[(1 + length(Z)) + 1]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muX <- eta0
    if (length(Z) > 0) {
      muX <- muX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
    }
    eX <- uncens_data[, X] - muX
    pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    # -------------------------------------- Calculate
  } else if (distX == "log-normal") {
    # Get parameters ---------------------------------
    eta_params <- params[-c(1:length(theta_params))]
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(Z))] }
    sigX <- eta_params[(1 + length(Z)) + 1]
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muX <- eta0
    if (length(Z) > 0) {
      muX <- muX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
    }
    pXgivZ <- (1 / (uncens_data[, X] * sigX * sqrt(2 * pi))) * exp(-1 * (log(uncens_data[, X]) - muX)^2 / (2 * sigX^2))
    # -------------------------------------- Calculate
  } else if (distX == "gamma") {
    # Get parameters ---------------------------------
    eta_params <- params[-c(1:length(theta_params))]
    eta0 <- eta_params[1]
    if (!is.null(Z)) { eta1 <- eta_params[2:(1 + length(Z))] }
    shapeX <- exp(eta_params[(1 + length(Z)) + 1])
    # --------------------------------- Get parameters
    # Calculate --------------------------------------
    muX <- eta0
    if (length(Z) > 0) {
      muX <- exp(muX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1)))
    }
    pXgivZ <- (1 / gamma(shapeX)) * ((shapeX / muX) ^ shapeX) * ((uncens_data[, X])^(shapeX - 1)) * exp(-1 * shapeX * uncens_data[, X] / muX)
    # -------------------------------------- Calculate
  } else if (distX == "inverse-gaussian") {
    # Get parameters ---------------------------------
    eta_params <- params[-c(1:length(theta_params))]
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
    eta_params <- params[-c(1:length(theta_params))]
    shapex <- eta_params[1]; eta_params <- eta_params[-1]
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
    eta_params <- params[-c(1:length(theta_params))]
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
    eta_params <- params[-c(1:length(theta_params))]
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
  
  ####################################################
  # Calculate joint density P(Y,X,Z) #################
  ####################################################
  uncens_data <- data.frame(cbind(uncens_data, jointP = pYgivXZ * pXgivZ))
  ####################################################
  # Calculate the log-likelihood #####################
  ####################################################
  # Log-likelihood contribution of uncensored X ------
  ll <- sum(log(uncens_data[, "jointP"]))
  # ------ Log-likelihood contribution of uncensored X
  # Log-likelihood contribution of censored X --------
  joint_dens <- function(x, Yi, Zi) {
    ####################################################
    # Analysis model P(Y|X,Z) ##########################
    ####################################################
    if (distY == "normal") {
      # Calculate --------------------------------------
      muY <- beta0 + beta1 * matrix(data = x, ncol = 1)
      if (length(Z) > 0) {
        muY <- muY + as.numeric(data.matrix(Zi) %*% matrix(data = beta2, ncol = 1))
      }
      muY <- data.matrix(muY)
      eY <- as.numeric(Yi) - muY
      pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
      # -------------------------------------- Calculate
    } else if (distY == "binomial") {
      # Calculate --------------------------------------
      muY <- beta0 + beta1 * matrix(data = x, ncol = 1)
      if (length(Z) > 0) {
        muY <- muY + as.numeric(data.matrix(Zi) %*% matrix(data = beta2, ncol = 1))
      }
      muY <- data.matrix(muY)
      pYgivXZ <- exp(- (1 - Yi) * muY) / (1 + exp(- muY))
      # -------------------------------------- Calculate
    }
    ####################################################
    # Predictor model P(X|Z) ###########################
    ####################################################
    if (distX == "normal") {
      # Calculate --------------------------------------
      muX <- eta0
      if (length(Z) > 0) {
        muX <- muX + as.numeric(data.matrix(Zi) %*% matrix(data = eta1, ncol = 1))
      }
      eX <- x - muX
      pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
      # -------------------------------------- Calculate
    } else if (distX == "log-normal") {
      # Calculate --------------------------------------
      muX <- eta0
      if (length(Z) > 0) {
        muX <- muX + as.numeric(data.matrix(Zi) %*% matrix(data = eta1, ncol = 1))
      }
      pXgivZ <- (1 / (x * sigX * sqrt(2 * pi))) * exp(-1 * (log(x) - muX)^2 / (2 * sigX^2))
      # -------------------------------------- Calculate
    } else if (distX == "gamma") {
      # Calculate --------------------------------------
      muX <- eta0
      if (length(Z) > 0) {
        muX <- exp(muX + as.numeric(data.matrix(Zi) %*% matrix(data = eta1, ncol = 1)))
      }
      pXgivZ <- (1 / gamma(shapeX)) * ((shapeX / muX) ^ shapeX) * ((x)^(shapeX - 1)) * exp(-1 * shapeX * x / muX)
      # -------------------------------------- Calculate
    } else if (distX == "inverse-gaussian") {
      # Calculate --------------------------------------
      muX <- eta0
      if (length(Z) > 0) {
        muX <- exp(muX + as.numeric(data.matrix(Zi) %*% matrix(data = eta1, ncol = 1)))
      }
      pXgivZ <- sqrt((shapeX / (2 * pi * x^3))) * exp(-1 * (shapeX * (x - muX)^2) / (2 * muX^2 * x))
      # -------------------------------------- Calculate
    } else if (distX == "weibull") {
      # Calculate --------------------------------------
      shapex <- eta0
      pXgivZ <- (shapex / scalex) * ((x / scalex) ^ (shapex - 1)) * exp(-1 * (x / scalex)^ shapex)
      # -------------------------------------- Calculate
    } else if (distX == "exponential") {
      # Calculate --------------------------------------
      lambdaX <- eta0
      if (length(Z) > 0) {
        lambdaX <- lambdaX + as.numeric(data.matrix(Zi) %*% matrix(data = eta1, ncol = 1))
      }
      pXgivZ <- lambdaX * exp(- lambdaX * x)
      # -------------------------------------- Calculate
    } else if (distX == "poisson") {
      muX <- eta0
      if (length(Z) > 0) {
        muX <- muX + as.numeric(data.matrix(Zi) %*% matrix(data = eta1, ncol = 1))
      }
      pXgivZ <- ((muX^x) * (exp(-1 * muX))) / factorial(x)
      # -------------------------------------- Calculate
    }
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
  ll <- ll + sum(log_integral)
  # -------- Log-likelihood contribution of censored X
  # Return (-1) x log-likelihood for use with nlm() --
  return(- ll)
  # -- Return (-1) x log-likelihood for use with nlm()
}
