#' Maximum likelihood estimator (MLE) for censored predictor in generalized linear models (GLM)
#'
#' @param Y name of outcome variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param D name of event indicator, defined to be \code{= 1} if \code{X} was uncensored and \code{0} otherwise.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data a dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but other options are \code{"bernoulli"}, \code{"lognormal"}, \code{"gamma"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"lognormal"}, \code{"gamma"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
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

  ## Create a list containing data and other inputs
  ### This list will be fed to functions later to reduce the number of arguments
  dataObj = list(Y = Y,
                 W = W,
                 D = D,
                 Z = Z,
                 data = data,
                 rightCens = rightCens, # direction of censoring (if TRUE, right; if FALSE, left)
                 distY = distY, # distribution for Y|X,Z
                 distX = distX # distribution for Z|Z
                 )

  # class(dataObj) = c(paste0(distY, "Y"), # distribution for Y|X,Z
  #                    paste0(distX, "X"), # distribution for Z|Z
  #                    ifelse(test = rightCens, # direction of censoring
  #                           yes = "rightCensRd",
  #                           no = "leftCensRd"
  #                           )
  #                    )

  # Initial parameter values
  params0 = init_vals_cc(dataObj = dataObj,
                         D = D,
                         steptol = steptol,
                         iterlim = iterlim,
                         verbose = verbose)

  if (verbose) {
    print("Fit full MLE:")
  }
  mod = suppressWarnings(
    nlm(f = loglik,
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

    ## Augment dataObj with parameters at convergence and SEs
    conv_dataObj = dataObj
    conv_dataObj$params = param_est
    conv_dataObj$vcov = param_vcov
    conv_dataObj$se = param_se

    if (robcov) {
      # Derivatives of the log-likelihood
      first_deriv = calc_deriv_loglik(dataObj = conv_dataObj,
                                      subdivisions = subdivisions)
      second_deriv = calc_deriv2_loglik(dataObj = conv_dataObj,
                                        subdivisions = subdivisions)

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
      A_inv = tryCatch(expr = solve(A),
                       error = function(c) matrix(data = NA,
                                                  nrow = length(param_est),
                                                  ncol = length(param_est))
      )
      ## Sandwich covariance
      n = nrow(data)
      param_rob_vcov = A_inv %*% B %*% t(A_inv)
      param_rob_se = sqrt(diag(param_rob_vcov)) / sqrt(n)
    } else {
      param_rob_vcov = matrix(data = NA,
                              nrow = length(param_est),
                              ncol = length(param_est))
      param_rob_se = rep(NA,
                         length(param_est))
    }

    ## Augment dataObj with robust vcov and SEs
    conv_dataObj$rob_vcov = param_rob_vcov
    conv_dataObj$rob_se = param_rob_se
  } else {
    param_est = param_se = param_rob_se = rep(NA,
                                              times = length(mod$estimate))
    param_vcov = param_rob_vcov = matrix(data = NA,
                                         nrow = length(param_est),
                                         ncol = length(param_est))

    ## Augment dataObj with parameters at convergence and SEs
    conv_dataObj = dataObj
    conv_dataObj$params = param_est
    conv_dataObj$vcov = param_vcov
    conv_dataObj$se = param_se
    conv_dataObj$rob_vcov = param_rob_vcov
    conv_dataObj$rob_se = param_rob_se
  }

  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  ####################################################
  modY = get_modY(object = conv_dataObj)

  ## Remove outcome model parameters/standard errors
  dim_beta = length(get_beta_params(object = conv_dataObj))
  conv_dataObj$params = conv_dataObj$params[-c(1:dim_beta)]
  conv_dataObj$se = conv_dataObj$se[-c(1:dim_beta)]
  conv_dataObj$rob_se = conv_dataObj$rob_se[-c(1:dim_beta)]

  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  modX = get_modX(object = conv_dataObj)

  ####################################################
  # Return model results #############################
  ####################################################
  ## Start with a list of class "glmCensRd"
  res = list(outcome_model = unclass(modY),
             covariate_model = unclass(modX),
             code = mod$code,
             vcov = param_vcov,
             rvcov = param_rob_vcov)
  class(res) = "glmCensRd"

  return(res)
}
