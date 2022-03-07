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

  if (distY == "normal") {
    params0 <- c(rep(0, length(c(1, X, Z))), var(data[, Y]))
  }

  if (distX %in% c("normal", "log-normal")) {
    params0 <- c(params0, rep(0, length(c(1, Z))), var(data[, X], na.rm = TRUE))
  } else if (distX %in% c('gamma', "inverse-gaussian", "weibull")) {
    params0 <- c(params0, 0.1, rep(0.1, length(c(1, Z))))
  }

  suppressWarnings(
    mod <- nlm(f = loglik, p = params0, steptol = steptol, iterlim = iterlim, hessian = TRUE,
               Y = Y, X = X, D = D, W = W, Z = Z, partX = partX, distY = distY, distX = distX, data = data)
  )
  param_est <- mod$estimate
  #param_se <- sqrt(diag(solve(mod$hessian)))

  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  ####################################################
  if (distY == "normal") {
    modY_est <- param_est[1:(length(c(X, Z)) + 1)]
    modY_se <- NA # param_se[1:(length(c(X, Z)) + 1)]
    modY_sigma2 <- param_est[(length(c(X, Z)) + 1) + 1] ^ 2
    modY_coeff <- data.frame(coeff = modY_est, se = modY_se)
    rownames(modY_coeff) <- c("(Intercept)", X, Z)
    modY <- list(distY = distY, mean = modY_coeff, sigma2 = modY_sigma2)
    param_est <- param_est[-c(1:(length(c(X, Z)) + 2))]
    #param_se <- param_se[-c(1:(length(c(X, Z)) + 2))]
  } else if (distY == "binomial") {
    modY_coeff <- param_est[1:(length(c(X, Z)) + 1)]
    modY_se <- NA #param_se[1:(length(c(X, Z)) + 1)]
    modY_coeff <- data.frame(coeff = modY_est, se = modY_se)
    rownames(modY_coeff) <- c("(Intercept)", X, Z)
    modY <- list(distY = distY, mean = modY_coeff)
    param_est <- param_est[-c(1:(length(c(X, Z)) + 1))]
    #param_se <- param_se[-c(1:(length(c(X, Z)) + 1))]
  }

  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  if (distX %in% c("normal", "log-normal")) {
    modX_est <- param_est[1:(length(Z) + 1)]
    modX_se <- NA # param_se[1:(length(Z) + 1)]
    modX_sigma2 <- param_est[(length(Z) + 1) + 1] ^ 2
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- c("(Intercept)", Z)
    modX <- list(distX = distX, mean = modX_coeff, sigma2 = modX_sigma2)
  } else if (distX %in% c("gamma", "inverse-gaussian")) {
    modX_shape_est <- param_est[1]
    modX_shape_se <- NA #param_se[1:(length(Z) + 1)]
    modX_shape <- data.frame(coeff = modX_shape_est, se = modX_shape_se)
    rownames(modX_shape) <- c("(Intercept)")
    param_est <- param_est[-1]

    modX_est <- param_est[1:(length(Z) + 1)]
    modX_se <- NA # param_se[1:(length(Z) + 1)]
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- c("(Intercept)", Z)

    modX <- list(distX = distX, mean = modX_coeff, shape = modX_shape)
  } else if (distX == "weibull") {
    modX_shape_est <- param_est[1]
    modX_shape_se <- NA #param_se[1:(length(Z) + 1)]
    modX_shape <- data.frame(coeff = modX_shape_est, se = modX_shape_se)
    rownames(modX_shape) <- c("(Intercept)")
    param_est <- param_est[-1]

    modX_est <- param_est[1:(length(Z) + 1)]
    modX_se <- NA # param_se[1:(length(Z) + 1)]
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- rownames(modX_shape) <- c("(Intercept)", Z)

    modX <- list(distX = distX, scale = modX_coeff, shape = modX_shape)
  } else if (distX %in% c("exponential", "poisson")) {
    modX_est <- param_est[1:(length(Z) + 1)]
    modX_se <- NA # param_se[1:(length(Z) + 1)]
    modX_coeff <- data.frame(coeff = modX_est, se = modX_se)
    rownames(modX_coeff) <- c("(Intercept)", Z)
    modX <- list(distX = distX, rate = modX_coeff)
  }

  return(list(outcome_model = modY, predictor_model = modX, code = mod$code))
}
