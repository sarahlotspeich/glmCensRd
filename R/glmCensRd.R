#' Maximum likelihood estimator (MLE) for censored predictor in generalized linear models (GLM)
#'
#' @param Y Name of column variable.
#' @param W Name of observed (i.e., censored) version of \code{X}.
#' @param D Name of event indicator, defined to be = 1 if \code{X} was uncensored.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{NULL}
#' @param partX Size of partition of unobserved \code{X} for censored subjects. Default is \code{50}.
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY Distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}.
#' @param distX Distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}.
#' @param steptol (Fed to \code{nlm()}) A positive scalar providing the minimum allowable relative step length. Default is \code{1e-4}.
#' @param iterlim (Fed to \code{nlm()}) A positive integer specifying the maximum number of iterations to be performed before the program is terminated. Default is \code{100}.
#'
#' @return A list with the following elements:
#' \item{outcome_model}{A list containing details of the fitted model for the outcome.}
#' \item{predictor_model}{A list containing details of the fitted model for the predictor.}
#' \item{code}{An integer indicating why the optimization process terminated. See \code{?nlm} for details on values.}
#'
#' @export
#' @importFrom maxLik maxLik
#'
glmCensRd <- function(Y, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data, steptol = 1e-4, iterlim = 100) {
  # Subset data to relevant, user-specified columns
  data <- data[, c(Y, W, D, Z)]
  # Create variable X = W
  data <- cbind(data, X = data[, W])
  ## Make it NA for censored (D = 0)
  data[data[, D] == 0, "X"] <- NA
  ## Define column X variable name
  X <- "X"

  # Initial parameter values
  if (distY == "normal") {
    params0 <- c(rep(0, length(c(1, X, Z))), var(data[, Y]))
  } else if (distY == "binomial") {
    params0 <- c(rep(0, length(c(1, X, Z))))
  }
  if (distX %in% c("normal", "log-normal")) {
    params0 <- c(params0, rep(0, length(c(1, Z))), var(data[, X], na.rm = TRUE))
  } else if (distX %in% c('gamma', "inverse-gaussian")) {
    params0 <- c(params0, 0.1, rep(0.1, length(c(1, Z))))
  } else if (distX == "weibull") {
    params0 <- c(params0, 1, rep(0.1, length(c(1, Z))))
  } else if (distX %in% c("exponential", "poisson")) {
    params0 <- c(params0, rep(0.1, length(c(1, Z))))
  }

  suppressWarnings(
    mod <- nlm(f = loglik, p = params0, steptol = steptol, iterlim = iterlim, hessian = TRUE,
               Y = Y, X = X, D = D, W = W, Z = Z, partX = partX, distY = distY, distX = distX, data = data)
  )
  param_est <- mod$estimate
  #param_cov <- solve(mod$hessian)

  # Derivatives of the log-likelihood
  first_deriv <- calc_deriv_loglik(mle = param_est, Y = Y, X = X, D = D, W = W, Z = Z,
                                   partX = partX, distY = distY, distX = distX, data = data)

  second_deriv <- calc_deriv2_loglik(mle = param_est, Y = Y, X = X, D = D, W = W, Z = Z,
                                     partX = partX, distY = distY, distX = distX, data = data)

  # Sandwich covariance estimator
  ## Sandwich meat
  rep_each <- first_deriv[, rep(x = 1:ncol(first_deriv), each = length(param_est))]
  rep_times <- first_deriv[, rep(x = 1:ncol(first_deriv), times = length(param_est))]
  entriesB <- colMeans(x = rep_each * rep_times)
  B <- matrix(data = entriesB,
              nrow = length(param_est),
              ncol = length(param_est),
              byrow = TRUE)

  ## Sandwich bread
  entriesA <- colMeans(x = second_deriv)
  A <- matrix(data = entriesA,
              nrow = length(param_est),
              ncol = length(param_est),
              byrow = TRUE)

  ## Sandwich covariance
  n <- nrow(data)
  param_cov <- (solve(A) %*% B %*% t(solve(A)))
  param_se <- sqrt(diag(param_cov)) / sqrt(n)

  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  ####################################################
  if (distY == "normal") {
    dim_beta <- length(c(X, Z)) + 1
    modY_est <- param_est[1:dim_beta]
    modY_se <- param_se[1:dim_beta]
    modY_sigma2 <- param_est[dim_beta + 1] ^ 2
    modY_coeff <- data.frame(coeff = modY_est, se = modY_se)
    rownames(modY_coeff) <- c("(Intercept)", X, Z)
    modY <- list(distY = distY, mean = modY_coeff, sigma2 = modY_sigma2)
    param_est <- param_est[-c(1:(dim_beta + 1))]
    param_se <- param_se[-c(1:(dim_beta + 1))]
  } else if (distY == "binomial") {
    dim_beta <- length(c(X, Z)) + 1
    modY_est <- param_est[1:dim_beta]
    modY_se <- param_se[1:dim_beta]
    modY_coeff <- data.frame(coeff = modY_est, se = modY_se)
    rownames(modY_coeff) <- c("(Intercept)", X, Z)
    modY <- list(distY = distY, mean = modY_coeff)
    param_est <- param_est[-c(1:dim_beta)]
    param_se <- param_se[-c(1:dim_beta)]
  }

  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  if (distX %in% c("normal", "log-normal")) {
    dim_eta <- length(Z) + 1
    modX_est <- param_est[1:dim_eta]
    modX_se <- param_se[1:dim_eta]
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- c("(Intercept)", Z)
    modX_sigma2 <- param_est[dim_eta + 1] ^ 2
    modX <- list(distX = distX, mean = modX_coeff, sigma2 = modX_sigma2)
  } else if (distX %in% c("gamma", "inverse-gaussian")) {
    modX_shape_est <- param_est[1]
    modX_shape_se <- param_se[1]
    modX_shape <- data.frame(coeff = modX_shape_est, se = modX_shape_se)
    rownames(modX_shape) <- c("(Intercept)")
    param_est <- param_est[-1]
    param_se <- param_se[-1]

    dim_eta <- length(Z) + 1
    modX_est <- param_est[1:dim_eta]
    modX_se <- param_se[1:dim_eta]
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- c("(Intercept)", Z)

    modX <- list(distX = distX, mean = modX_coeff, shape = modX_shape)
  } else if (distX == "weibull") {
    modX_shape_est <- param_est[1]
    modX_shape_se <- param_se[1]
    modX_shape <- data.frame(coeff = modX_shape_est, se = modX_shape_se)
    rownames(modX_shape) <- c("(Intercept)")
    param_est <- param_est[-1]
    param_se <- param_se[-1]

    dim_eta <- (length(Z) + 1)
    modX_est <- param_est[1:dim_eta]
    modX_se <- param_se[1:dim_eta]
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- c("(Intercept)", Z)

    modX <- list(distX = distX, scale = modX_coeff, shape = modX_shape)
  } else if (distX %in% c("exponential", "poisson")) {
    dim_eta <- length(Z) + 1
    modX_est <- param_est[1:dim_eta]
    modX_se <- param_se[1:dim_eta]
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- c("(Intercept)", Z)
    modX <- list(distX = distX, rate = modX_coeff)
  }

  return(list(outcome_model = modY, predictor_model = modX, code = mod$code))
}
