#' Observed-data log-likelihood for normal linear regression with a censored predictor modeled with normal linear regression, also.

#' @param params0 Parameter values.
#' @param Y Name of outcome variable.
#' @param X Name of censored predictor variable.
#' @param W Name of observed (i.e., censored) version of \code{X}.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{NULL}
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{W}, \code{Z}.
#' @param steptol (Fed to \code{nlm()}) A positive scalar providing the minimum allowable relative step length. Default is \code{steptol = 1e-6}.
#' @param iterlim (Fed to \code{nlm()}) A positive integer specifying the maximum number of iterations to be performed before the program is terminated. Default is \code{iterlim = 100}.

#' @return A vector containing the approximate integral over the joint density of censored subjects from \code{complete_data_cens}.
#'
loglik_normal_normal <- function(params, Y, X, Z = NULL, complete_data) {
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
    muY <- beta0 + beta1 * complete_data[, X] + data.matrix(complete_data[, Z]) %*% beta2
  } else {
    muY <- beta0 + beta1 * complete_data[, X]
  }
  eY <- complete_data[, Y] - muY
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
    muX <- eta0 + data.matrix(complete_data[, Z]) %*% eta1
  } else {
    muX <- eta0
  }
  eX <- complete_data[, X] - muX
  pXgivZ <- 1 / sqrt(2 * pi * sigX ^ 2) * exp(- eX ^ 2 / (2 * sigX ^ 2))
  # ---------------------------------------- Calculate
  ####################################################
  # Calculate joint density P(Y,X,Z) #################
  ####################################################
  complete_data <- cbind(complete_data, jointP = pYgivXZ * pXgivZ)
  ####################################################
  # Calculate the log-likelihood #####################
  ####################################################
  # Log-likelihood contribution of uncensored X ------
  ll <- sum(log(complete_data[1:n1, "jointP"]))
  # ------ Log-likelihood contribution of uncensored X
  # Log-likelihood contribution of censored X --------
  log_integral <- log(integrateCensRd(X = X, complete_data_cens = complete_data[-c(1:n1), ]))
  log_integral[log_integral == -Inf] <- 0
  ll <- ll + sum(log_integral)
  # -------- Log-likelihood contribution of censored X
  # Return - log-likelihood for use with nlm() -------
  return(- ll)
}
