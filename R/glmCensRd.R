#' Maximum likelihood estimator (MLE) for generalized linear models (GLMs) with a censored covariate
#'
#' @param Y name of outcome variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param D name of event indicator, defined to be \code{= 1} if \code{X} was uncensored and \code{0} otherwise.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but other options are \code{"binomial"}, \code{"log-normal"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param cens type of censoring assumed for \code{X}. Default is \code{"right"}, but the other option is \code{"left"}.
#' @param robcov logical. If \code{TRUE} (the default), the robust sandwich covariance estimator of parameter standard errors is included in the output.
#' @param init_vals (optional) initial values for \code{outcome_model} and \code{covariate_model}. Default is \code{NULL} and naive initial values will be used. 
#' @param subdivisions (fed to \code{integrate}) the maximum number of subintervals used to integrate over unobserved \code{X} for censored subjects. Default is \code{100}.
#' @param steptol (fed to \code{nlm}) a positive scalar providing the minimum allowable relative step length. Default is \code{1E-6}.
#' @param iterlim (fed to \code{nlm}) a positive integer specifying the maximum number of iterations to be performed before the program is terminated. Default is \code{100}.
#'
#' @return A list with the following elements:
#' \item{outcome_model}{a list containing details of the fitted model for the outcome.}
#' \item{covariate_model}{a list containing details of the fitted model for the predictor.}
#' \item{code}{an integer indicating why the optimization process terminated. See \code{?nlm} for details on values.}
#' \item{vcov}{variance-covariance matrix of the model parameters.}
#' \item{rvcov}{robust sandwich variance-covariance matrix of the model parameters.}
#'
#' @export
#'
glmCensRd = function(Y, W, D, Z = NULL, data,  distY = "normal", distX = "normal", cens = "right", 
                     robcov = TRUE, init_vals = NULL, subdivisions = 100, steptol = 1E-6, iterlim = 100) {
  # Subset data to relevant, user-specified columns
  data = data[, c(Y, W, D, Z)]

  # Create variable X = W
  data = cbind(data, X = data[, W])

  ## Make it NA for censored (D = 0)
  data[data[, D] == 0, "X"] = NA

  ## Define column X variable name
  X = "X"

  # Initial parameter values
  # If not supplied by user, use defaults 
  if (is.null(init_vals)) {
    ## Naive initial parameters
    params0_n = init_vals(Z = Z,
                          distY = distY,
                          distX = distX)
    
    ## Use complete-case MLE
    cc_data = data[data[, D] == 1, ]
    cc_mod = nlm(f = loglik_cc,
                 p = params0_n,
                 Y = Y,
                 X = X,
                 Z = Z,
                 data = cc_data,
                 distY = distY,
                 distX = distX,
                 steptol = steptol,
                 iterlim = iterlim,
                 hessian = FALSE
    )
    
    ## Check for convergence 
    if (cc_mod$code <= 2 & cc_mod$iterations > 1) {
      init_vals = cc_mod$estimate ### If TRUE, return estimates
    } else {
      init_vals = params0_n ### If FALSE, return naive estimates instead
    }
  } 
  
  theta_bounds = bounds(Z = Z, 
                        distY = distY, 
                        distX = distX)
  
  mod = optim(par = init_vals, 
              fn = loglik, #gr = gloglik,
              Y = Y, 
              X = X,
              D = D, 
              W = W, 
              Z = Z, 
              data = data,
              distY = distY, 
              distX = distX, 
              cens = cens, 
              subdivisions = subdivisions,
              method = "L-BFGS-B",
              control = list(maxit = iterlim),
              lower = theta_bounds$lower,
              upper = theta_bounds$upper,
              hessian = T)

  # Check that optim() converged 
  if (mod$convergence == 0) {
    ## Parameter estimates
    est = mod$par
    numpar = length(est)
    
    ## Variance/standard error estimates
    vcov = tryCatch(expr = solve(mod$hessian),
                    error = function(c) 
                      matrix(data = NA,
                             nrow = numpar,
                             ncol = numpar)
    )
    se = sqrt(diag(vcov))
    
    ## Robust variance/standard error estimates
    if (robcov) {
      # Derivatives of the log-likelihood
      first_deriv = calc_deriv_loglik(params = est,
                                      Y = Y,
                                      X = X,
                                      D = D,
                                      W = W,
                                      Z = Z,
                                      subdivisions = subdivisions,
                                      distY = distY,
                                      distX = distX,
                                      data = data)
      
      second_deriv = calc_deriv2_loglik(params = est,
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
      rep_each = first_deriv[, rep(x = 1:ncol(first_deriv), each = numpar)]
      rep_times = first_deriv[, rep(x = 1:ncol(first_deriv), times = numpar)]
      entriesB = colMeans(x = rep_each * rep_times)
      B = matrix(data = entriesB,
                 nrow = numpar,
                 ncol = numpar,
                 byrow = TRUE)
      
      ## Sandwich bread
      entriesA = colMeans(x = second_deriv)
      A = matrix(data = entriesA,
                 nrow = numpar,
                 ncol = numpar,
                 byrow = TRUE)
      
      ## Build sandwich
      n = nrow(data)
      rob_vcov = solve(A) %*% B %*% t(solve(A)) / n
      rob_se = sqrt(diag(rob_vcov))
    } else {
      rob_vcov = matrix(data = NA,
                        nrow = numpar,
                        ncol = numpar)
      rob_se = rep(NA, numpar)
    }
  } else {
    est = se = rob_se = rep(NA, times = length(mod$estimate))
    vcov = rob_vcov = matrix(data = NA,
                             nrow = numpar,
                             ncol = numpar)
  }

  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  ####################################################
  summ = summarize_models(est = est, 
                          se = se, 
                          rob_se = rob_se, 
                          X = X,
                          Z = Z, 
                          distY = distY, 
                          distX = distX)

  # Return model results in a list
  return(list(outcome_model = summ$modY,
              covariate_model = summ$modX,
              code = mod$convergence,
              vcov = vcov,
              rvcov = rob_vcov
              )
         )
}
