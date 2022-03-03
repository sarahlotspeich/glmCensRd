#' Maximum likelihood estimator (MLE) for censored predictor in generalized linear models (GLM)
#'
#' @param params0 Starting parameter values.
#' @param Y Name of column variable.
#' @param X Name of censored covariate variable.
#' @param W Name of observed (i.e., censored) version of \code{X}.
#' @param D Name of event indicator, defined to be = 1 if \code{X} was uncensored.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{NULL}
#' @param partX Size of partition of unobserved \code{X} for censored subjects. Default is \code{50}.
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY Distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}.
#' @param distX Distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}.
#' @param steptol (Fed to \code{nlm()}) A positive scalar providing the minimum allowable relative step length. Default is \code{1e-6}.
#' @param iterlim (Fed to \code{nlm()}) A positive integer specifying the maximum number of iterations to be performed before the program is terminated. Default is \code{100}.
#'
#' @return A list with the following two elements:
#' \item{coeff}{A dataframe containing coefficient and standard error estimates at convergence.}
#' \item{code}{An integer indicating why the optimization process terminated. See \code{?nlm} for details on values.}
#'
#' @export
#'
glmCensRd <- function(params0, Y, X, W, D, Z = NULL, partX = 50, distY = "normal", distX = "normal", data, steptol = 1e-6, iterlim = 100) {
  suppressWarnings(
    mod <- nlm(f = loglik, p = params0, steptol = steptol, iterlim = iterlim, hessian = TRUE,
               Y = Y, X = X, D = D, W = W, Z = Z, partX = partX, distY = distY, distX = distX, data = data)
  )
  #return(mod)
  param_est <- mod$estimate
  param_se <- sqrt(diag(solve(mod$hessian)))
  param_df <- data.frame(est = param_est, se = param_se)
  rownames(param_df) <- c(paste0("beta", 0:length(c(X, Z))), "sigmaY", paste0("eta", 0:length(c(X, Z))), "sigmaX")
  return(list(coeff = param_df, code = mod$code))
}
