#' Maximum likelihood estimator (MLE) for censored predictor in generalized linear models (GLM)
#'
#' @param Y name of outcome variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param D name of event indicator, defined to be \code{= 1} if \code{X} was uncensored and \code{0} otherwise.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data a dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but other options are \code{"bernoulli"}, \code{"log-normal"}, \code{"gamma"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param rightCens logical. If \code{TRUE} (the default), right censoring is assumed. If \code{FALSE}, left censoring is assumed.
#' @param robcov logical. If \code{TRUE} (the default), the robust sandwich covariance estimator of parameter standard errors is included in the output.
#' @param subdivisions (fed to \code{integrate}) the maximum number of subintervals used to integrate over unobserved \code{X} for censored subjects. Default is \code{100}.
#' @param steptol (fed to \code{nlm}) a positive scalar providing the minimum allowable relative step length. Default is \code{1E-6}.
#' @param iterlim (fed to \code{nlm}) a positive integer specifying the maximum number of iterations to be performed before the program is terminated. Default is \code{100}.
#' @param verbose logical. If \code{TRUE}, text is printed to show progress. Default is \code{FALSE}.
#'
#' @return A list with the following elements:
#' \item{outcome_model}{a list containing details of the fitted model for the outcome.}
#' \item{predictor_model}{a list containing details of the fitted model for the predictor.}
#' \item{code}{an integer indicating why the optimization process terminated. See \code{?nlm} for details on values.}
#' \item{vcov}{variance-covariance matrix of the model parameters.}
#' \item{rvcov}{robust sandwich variance-covariance matrix of the model parameters.}
#'
#' @export
#'
glmCensRd = function(Y, W, D, Z = NULL, data,  distY = "normal", distX = "normal",
                     rightCens = TRUE, robcov = TRUE, subdivisions = 100, steptol = 1E-6,
                     iterlim = 100, verbose = FALSE) {
  # Subset data to relevant, user-specified columns
  data = data[, c(Y, W, D, Z)]

  # Create variable X = W
  data = cbind(data, X = data[, W])

  ## Make it NA for censored (D = 0)
  data[data[, D] == 0, "X"] = NA

  ## Define column X variable name
  X = "X"

  ## Create a list containing data and other inputs
  ### This list will be fed to functions later to reduce the number of arguments
  dataObj = list(Y = Y,
                 W = W,
                 D = D,
                 Z = Z,
                 data = data,
                 rightCens = rightCens
                 )

  class(dataObj) = c("dataCensRd", # data for glmCensRd
                     paste0(distY, "Y"), # distribution for Y|X,Z
                     paste0(distX, "X"), # distribution for Z|Z
                     ifelse(test = rightCens, # direction of censoring
                            yes = "rightCensRd",
                            no = "leftCensRd"
                            )
                     )

  # Initial parameter values
  ## Naive initial params (for complete-case)
  if (verbose) {
    print("Start with naive initial values:")
  }
  params0 = init_vals(object = dataObj)
  if (verbose) {
    print(params0)
  }
  if (verbose) {
    print("Fit complete-case initial values:")
  }

  ## Use complete-case MLE
  cc_dataObj = dataObj
  cc_dataObj$data = data[data[, D] == 1, ]
  suppressWarnings(
    cc_mod <- nlm(f = cc_loglik,
                  p = params0,
                  dataObj = cc_dataObj,
                  steptol = steptol,
                  iterlim = iterlim,
                  hessian = FALSE
                  )
    )
  if (cc_mod$code <= 2 & cc_mod$iterations > 1) {
    params0 = cc_mod$estimate
  }
  if (verbose) {
    print(params0)
  }

  if (verbose) {
    print("Fit full MLE:")
  }
  suppressWarnings(
    mod <- nlm(f = loglik,
               p = params0,
               dataObj = dataObj,
               subdivisions = subdivisions,
               steptol = steptol,
               iterlim = iterlim,
               hessian = TRUE)
  )
  if (verbose) {
    print(mod$estimate)
  }

  # Check that nlm() actually iterated
  if (mod$iterations > 1) {
    param_est = mod$estimate
    param_vcov = tryCatch(expr = solve(mod$hessian),
                           error = function(c) matrix(data = NA,
                                                      nrow = length(param_est),
                                                      ncol = length(param_est))
    )
    param_se = sqrt(diag(param_vcov))

    if (robcov) {
      # Derivatives of the log-likelihood
      first_deriv = calc_deriv_loglik(params = param_est,
                                       Y = Y,
                                       X = X,
                                       D = D,
                                       W = W,
                                       Z = Z,
                                       subdivisions = subdivisions,
                                       distY = distY,
                                       distX = distX,
                                       data = data)

      second_deriv = calc_deriv2_loglik(params = param_est,
                                         Y = Y,
                                         X = X,
                                         D = D,
                                         W = W,
                                         Z = Z,
                                         subdivisions = subdivisions,
                                         distY = distY,
                                         distX = distX,
                                         data = data)

      # Sandwich covariance estimator
      ## Sandwich meat
      rep_each = first_deriv[, rep(x = 1:ncol(first_deriv), each = length(param_est))]
      rep_times = first_deriv[, rep(x = 1:ncol(first_deriv), times = length(param_est))]
      entriesB = colMeans(x = rep_each * rep_times)
      B = matrix(data = entriesB,
                  nrow = length(param_est),
                  ncol = length(param_est),
                  byrow = TRUE)

      ## Sandwich bread
      entriesA = colMeans(x = second_deriv)
      A = matrix(data = entriesA,
                  nrow = length(param_est),
                  ncol = length(param_est),
                  byrow = TRUE)

      ## Sandwich covariance
      n = nrow(data)
      param_rob_vcov = solve(A) %*% B %*% t(solve(A))
      param_rob_se = sqrt(diag(param_rob_vcov)) / sqrt(n)
    } else {
      param_rob_vcov = matrix(data = NA,
                               nrow = length(param_est),
                               ncol = length(param_est))
      param_rob_se = rep(NA, length(param_est))
    }
  } else {
    param_est = param_se = param_rob_se = rep(NA, times = length(mod$estimate))
    param_vcov = param_rob_vcov = matrix(data = NA,
                                           nrow = length(param_est),
                                           ncol = length(param_est))
  }

  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  ####################################################
  if (distY == "normal") {
    # Create coefficients dataframe for outcome model
    ## Mean parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_mean_est = param_est[1:dim_beta]
    modY_mean_se = param_se[1:dim_beta]
    modY_mean_rse = param_rob_se[1:dim_beta]
    modY_mean = data.frame(coeff = modY_mean_est,
                            se = modY_mean_se,
                            robse = modY_mean_rse)
    rownames(modY_mean) = c("(Intercept)", X, Z)

    ## Error variance parameter (estimated directly)
    modY_sigma2 = param_est[dim_beta + 1] ^ 2

    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                 mean = modY_mean,
                 sigma2 = modY_sigma2)
  } else if (distY == "bernoulli") {
    # Create coefficients dataframe for outcome model
    ## Mean parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_mean_est = param_est[1:dim_beta]
    modY_mean_se = param_se[1:dim_beta]
    modY_mean_rse = param_rob_se[1:dim_beta]
    modY_mean = data.frame(coeff = modY_mean_est,
                            se = modY_mean_se,
                            robse = modY_mean_rse)
    rownames(modY_mean) = c("(Intercept)", X, Z)

    # Construct contents of "predictor_model" slot
    modY = list(distY = distY,
                 mean = modY_mean)
  } else if (distY %in% c("gamma", "inverse-gaussian")) {
    # Create coefficients dataframe for predictor model
    ## Shape parameter (estimated directly)
    modY_shape_est = param_est[1]
    modY_shape_se = param_se[1]
    modY_shape_rse = param_rob_se[1]
    modY_shape = data.frame(coeff = modY_shape_est,
                             se = modY_shape_se,
                             robse = modY_shape_rse)
    rownames(modY_shape) = c("(Intercept)")
    param_est = param_est[-1]
    param_se = param_se[-1]
    param_rob_se = param_rob_se[-1]

    ## Mean parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_mean_est = param_est[1:dim_beta]
    modY_mean_se = param_se[1:dim_beta]
    modY_mean_rse = param_rob_se[1:dim_beta]
    modY_mean = data.frame(coeff = modY_mean_est,
                            se = modY_mean_se,
                            robse = modY_mean_rse)
    rownames(modY_mean) = c("(Intercept)", X, Z)

    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                 mean = modY_mean,
                 shape = modY_shape)
  } else if (distY == "weibull") {
    ## Shape parameter (estimated directly)
    modY_shape_est = param_est[1]
    modY_shape_se = param_se[1]
    modY_shape_rse = param_rob_se[1]
    modY_shape = data.frame(coeff = modY_shape_est,
                             se = modY_shape_se,
                             robse = modY_shape_rse)
    rownames(modY_shape) = c("(Intercept)")
    param_est = param_est[-1]
    param_se = param_se[-1]
    param_rob_se = param_rob_se[-1]

    ## Scale parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_scale_est = param_est[1:dim_beta]
    modY_scale_se = param_se[1:dim_beta]
    modY_scale_rse = param_rob_se[1:dim_beta]
    modY_scale = data.frame(coeff = modY_scale_est,
                             se = modY_scale_se,
                             robse = modY_scale_rse)
    rownames(modY_scale) = c("(Intercept)", X, Z)

    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                 scale = modY_scale,
                 shape = modY_shape)
  } else if (distY %in% c("exponential", "poisson")) {
    # Create coefficients dataframe for outcome model
    ## Rate parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_rate_est = param_est[1:dim_beta]
    modY_rate_se = param_se[1:dim_beta]
    modY_rate_rse = param_rob_se[1:dim_beta]
    modY_rate = data.frame(coeff = modY_rate_est,
                            se = modY_rate_se,
                            robse = modY_rate_rse)
    rownames(modY_rate) = c("(Intercept)", X, Z)

    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                 rate = modY_rate)
  }
  # Remove outcome model parameters/standard errors
  param_est = param_est[-c(1:(dim_beta + 1))]
  param_se = param_se[-c(1:(dim_beta + 1))]
  param_rob_se = param_rob_se[-c(1:(dim_beta + 1))]

  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  if (distX %in% c("normal", "log-normal")) {
    # Create coefficients dataframe for predictor model
    ## Mean parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_mean_est = param_est[1:dim_eta]
    modX_mean_se = param_se[1:dim_eta]
    modX_mean_rse = param_rob_se[1:dim_eta]
    modX_mean = data.frame(coeff = modX_mean_est,
                            se = modX_mean_se,
                            robse = modX_mean_rse)
    rownames(modX_mean) = c("(Intercept)", Z)

    ## Error variance parameter (estimated directly)
    modX_sigma2 = param_est[dim_eta + 1] ^ 2

    # Construct contents of "predictor_model" slot
    modX = list(distX = distX,
                 mean = modX_mean,
                 sigma2 = modX_sigma2)
  } else if (distX %in% c("gamma", "inverse-gaussian")) {
    # Create coefficients dataframe for predictor model
    ## Shape parameter (estimated directly)
    modX_shape_est = param_est[1]
    modX_shape_se = param_se[1]
    modX_shape_rse = param_rob_se[1]
    modX_shape = data.frame(coeff = modX_shape_est,
                             se = modX_shape_se,
                             robse = modX_shape_rse)
    rownames(modX_shape) = c("(Intercept)")
    param_est = param_est[-1]
    param_se = param_se[-1]
    param_rob_se = param_rob_se[-1]

    ## Mean parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_mean_est = param_est[1:dim_eta]
    modX_mean_se = param_se[1:dim_eta]
    modX_mean_rse = param_rob_se[1:dim_eta]
    modX_mean = data.frame(coeff = modX_mean_est,
                            se = modX_mean_se,
                            robse = modX_mean_rse)
    rownames(modX_mean) = c("(Intercept)", Z)

    # Construct contents of "predictor_model" slot
    modX = list(distX = distX,
                 mean = modX_mean,
                 shape = modX_shape)
  } else if (distX == "weibull") {
    ## Shape parameter (estimated directly)
    modX_shape_est = param_est[1]
    modX_shape_se = param_se[1]
    modX_shape_rse = param_rob_se[1]
    modX_shape = data.frame(coeff = modX_shape_est,
                             se = modX_shape_se,
                             robse = modX_shape_rse)
    rownames(modX_shape) = c("(Intercept)")
    param_est = param_est[-1]
    param_se = param_se[-1]
    param_rob_se = param_rob_se[-1]

    ## Scale parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_scale_est = param_est[1:dim_eta]
    modX_scale_se = param_se[1:dim_eta]
    modX_scale_rse = param_rob_se[1:dim_eta]
    modX_scale = data.frame(coeff = modX_scale_est,
                             se = modX_scale_se,
                             robse = modX_scale_rse)
    rownames(modX_scale) = c("(Intercept)", Z)

    # Construct contents of "predictor_model" slot
    modX = list(distX = distX,
                 scale = modX_scale,
                 shape = modX_shape)
  } else if (distX %in% c("exponential", "poisson")) {
    ## Rate parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_rate_est = param_est[1:dim_eta]
    modX_rate_se = param_se[1:dim_eta]
    modX_rate_rse = param_rob_se[1:dim_eta]
    modX_rate = data.frame(coeff = modX_rate_est,
                            se = modX_rate_se,
                            robse = modX_rate_rse)
    rownames(modX_rate) = c("(Intercept)", Z)

    # Construct contents of "predictor_model" slot
    modX = list(distX = distX,
                 rate = modX_rate)
  }

  # Return model results in a list
  return(list(outcome_model = modY,
              predictor_model = modX,
              code = mod$code,
              vcov = param_vcov,
              rvcov = param_rob_vcov
              )
         )
}
