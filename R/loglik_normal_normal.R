#' Observed-data log-likelihood for normal linear regression with a censored predictor modeled with normal linear regression, also.

#' @param params0 Parameter values.
#' @param Y Name of outcome variable.
#' @param X Name of censored predictor variable.
#' @param W Name of observed (i.e., censored) version of \code{X}.
#' @param D Name of event indicator, defined to be = 1 if \code{X} was uncensored.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{NULL}
#' @param partX Size of partition of unobserved \code{X} for censored subjects. Default is \code{50}.
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{W}, \code{D}, and \code{Z}.

#' @return A vector containing the approximate integral over the joint density of censored subjects from \code{complete_data_cens}.
#'
loglik_normal_normal <- function(params, Y, X, W, D, Z = NULL, partX = 50, data) {
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
  # Get parameters -----------------------------------
  beta0 <- params[1]; params <- params[-1]
  beta1 <- params[1]; params <- params[-1]
  beta2 <- matrix(params[1:length(Z)], ncol = 1); params <- params[-c(1:length(Z))]
  sigY <- params[1]; params <- params[-1]
  # ----------------------------------- Get parameters
  # Calculate ----------------------------------------
  if (length(Z) > 1) {
    muY <- beta0 + beta1 * uncens_data[, X] + data.matrix(uncens_data[, Z]) %*% beta2
  } else {
    muY <- beta0 + beta1 * uncens_data[, X]
  }
  eY <- uncens_data[, Y] - muY
  pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
  # ---------------------------------------- Calculate
  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  # Get parameters -----------------------------------
  eta0 <- params[1]; params <- params[-1]
  if (!is.null(Z)) { eta1 <- params[1:length(Z)]; params <- params[-c(1:length(Z))] }
  sigX <- params[1]; params <- params[-1]
  # ----------------------------------- Get parameters
  # Calculate ----------------------------------------
  if (length(Z) > 1) {
    muX <- eta0 + data.matrix(uncens_data[, Z]) %*% eta1
  } else {
    muX <- eta0
  }
  eX <- uncens_data[, X] - muX
  pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
  # ---------------------------------------- Calculate
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
    # Calculate ----------------------------------------
    if (length(Z) > 1) {
      muY <- beta0 + beta1 * matrix(data = x, ncol = 1) + data.matrix(Zi) %*% beta2
    } else {
      muY <- beta0 + beta1 * matrix(data = x, ncol = 1)
    }
    muY <- data.matrix(muY)
    eY <- as.numeric(Yi) - muY
    pYgivXZ <- 1 / sqrt(2 * pi * sigY ^ 2) * exp(- eY ^ 2 / (2 * sigY ^ 2))
    # ---------------------------------------- Calculate
    ####################################################
    # Predictor model P(X|Z) ###########################
    ####################################################
    # Calculate ----------------------------------------
    if (length(Z) > 1) {
      muX <- eta0 + data.matrix(Zi) %*% eta1
    } else {
      muX <- eta0
    }
    eX <- x - muX
    pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
    return(pYgivXZ * pXgivZ)
  }
  integrate_joint_dens <- function(data_row) {
    data_row <- data.frame(t(data_row))
    joint_dens(x = seq(data_row[, W], data_row[, W] + 1, by = 0.25), Yi = data_row[Y], Zi = data_row[, Z])
    return(
      integrate(f = joint_dens, lower = data_row[, W], upper = Inf, subdivisions = partX,
              Yi = data_row[Y], Zi = data_row[, Z])$value
    )
  }
  integral <- apply(X = cens_data, MARGIN = 1, FUN = integrate_joint_dens)
  log_integral <- log(integral)
  log_integral[log_integral == -Inf] <- 0
  ll <- ll + sum(log_integral)
  # -------- Log-likelihood contribution of censored X
  # Return - log-likelihood for use with nlm() -------
  return(- ll)
}
