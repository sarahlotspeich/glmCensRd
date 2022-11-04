#' Observed-data log-likelihood

#' @param params parameter values.
#' @param Y name of outcome variable.
#' @param X name of censored predictor variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param D name of event indicator, defined to be \code{= 1} if \code{X} was uncensored and \code{0} otherwise.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param subdivisions (fed to \code{integrate}) the maximum number of subintervals used to integrate over unobserved \code{X} for censored subjects. Default is \code{100}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but \code{"binomial"} is the other option.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"inverse-gaussian"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#' @param cens type of censoring assumed for \code{X}. Default is \code{"right"}, but the other option is \code{"left"}.
#'
#' @export
#'
#' @return A scalar of the log-likelihood function (negated for use with \code{nlm}()).
#'
loglik <- function(params, Y, X, W, D, Z = NULL, data, subdivisions = 100, distY = "normal", distX = "normal", cens = "right") {
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
  # Joint density P(Y,X,Z) ###########################
  ####################################################
  pYXandZ_uncens <- calc_pYXandZ(x = uncens_data[, X],
                                 y = uncens_data[, Y],
                                 z = uncens_data[, Z],
                                 lengthZ = length(Z),
                                 distY = distY,
                                 distX = distX,
                                 params = params)

  ####################################################
  # Likelihood (Uncensored) ##########################
  ####################################################
  if (any(is.na(pYXandZ_uncens))) {
    # If params are out of domain, calc_pYXandZ returns NA
    ## And the log-likelihood needs to be arbitrarily "huge"
    return(1E8)
  } else {
    # Replace P(Y,X,Z) = 0 with P(Y,X,Z) = 1 so that
    ## log P(Y,X,Z) = 0.
    pYXandZ_uncens[pYXandZ_uncens == 0] = 1
    ll <- sum(log(pYXandZ_uncens))
  }

  ####################################################
  # Likelihood (Censored) ############################
  ####################################################
  if (nrow(cens_data) > 0) {
    integrate_pYXandZ <- function(data_row) {
      Wi <- as.numeric(data_row[W])
      Yi <- data_row[Y]
      Zi <- data_row[Z]
      return(
        tryCatch(expr = integrate(f = calc_pYXandZ,
                                  lower = ifelse(test = cens == "right", Wi, -Inf),
                                  upper = ifelse(test = cens == "right", Inf, Wi),
                                  subdivisions = subdivisions,
                                  y = Yi,
                                  z = Zi,
                                  lengthZ = length(Z),
                                  distY = distY,
                                  distX = distX,
                                  params = params)$value,
                 error = function(err) {0})
      )
    }
    int_pYXandZ_cens <- apply(X = cens_data,
                              MARGIN = 1,
                              FUN = integrate_pYXandZ)
    log_int_pYXandZ_cens <- log(int_pYXandZ_cens)
    log_int_pYXandZ_cens[log_int_pYXandZ_cens == -Inf] <- 0
    ll <- ll + sum(log_int_pYXandZ_cens)
  }

  # Return (-1) x log-likelihood for use with nlm() --
  return(- ll)
  # -- Return (-1) x log-likelihood for use with nlm()
}
