#' Gradient of observed-data log-likelihood for GLMs with a censored covariate

#' @param params parameter values.
#' @param Y name of outcome variable.
#' @param X name of censored predictor variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param D name of event indicator, defined to be \code{= 1} if \code{X} was uncensored and \code{0} otherwise.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but \code{"binomial"} is the other option.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"inverse-gaussian"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param cens type of censoring assumed for \code{X}. Default is \code{"right"}, but the other option is \code{"left"}.
#' @param subdivisions (fed to \code{integrate}) the maximum number of subintervals used to integrate over unobserved \code{X} for censored subjects. Default is \code{100}.
#'
#' @export
#'
#' @return A vector of the gradient of the log-likelihood function with respect to each parameter
#'
gloglik = function(params, Y, X, W, D, Z = NULL, data, distY = "normal", distX = "normal", cens = "right", subdivisions = 100) {
  
}
