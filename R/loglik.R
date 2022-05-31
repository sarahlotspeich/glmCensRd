#' Observed-data log-likelihood

#' @param params parameter values.
#' @param Y name of outcome variable.
#' @param X name of censored predictor variable.
#' @param W name of observed (censored) version of \code{X}.
#' @param D name of event indicator, defined to be \code{= 1} if \code{X} was uncensored and \code{0} otherwise.
#' @param Z (optional) name(s) of additional fully observed covariates. If none, \code{Z = NULL} (the default).
#' @param data a dataframe containing at least columns \code{Y}, \code{X}, \code{C}, \code{Z}.
#' @param subdivisions (fed to \code{integrate}) the maximum number of subintervals used to integrate over unobserved \code{X} for censored subjects. Default is \code{100}.
#' @param distY distribution assumed for \code{Y} given \code{X} and \code{Z}. Default is \code{"normal"}, but \code{"binomial"} is the other option.
#' @param distX distribution assumed for \code{X} given \code{Z}. Default is \code{"normal"}, but other options are \code{"log-normal"}, \code{"gamma"}, \code{"inverse-gaussian"}, \code{"weibull"}, \code{"exponential"}, or \code{"poisson"}.
#'
#' @export
#'
#' @return A scalar of the log-likelihood function (negated for use with \code{nlm}()).
#'
loglik <- function(params, Y, X, W, D, Z = NULL, data, subdivisions = 100, distY = "normal", distX = "normal") {
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
  if (distY == "normal") {
    # Subset parameters ------------------------------
    beta_params <- params[1:(3 + length(Z))]
    # ------------------------------ Subset parameters
  } else if (distY == "binomial") {
    # Subset parameters ------------------------------
    beta_params <- params[1:(2 + length(Z))]
    # ------------------------------ Subset parameters
  }
  pYgivXZ <- calc_pYgivXandZ(y = uncens_data[, Y],
                             x = uncens_data[, X],
                             z = uncens_data[, Z],
                             distY = distY,
                             beta_params = beta_params)

  ####################################################
  # Predictor model P(X|Z) ###########################
  ####################################################
  # Subset parameters --------------------------------
  eta_params <- params[-c(1:length(beta_params))]
  # -------------------------------- Subset parameters
  # Check for parameters outside domain --------------
  if (distX == "gamma") {
    # Shape and scale of Gamma > 0 -------------------
    shapeX <- eta_params[1]
    meanX <- eta_params[2]
    if (length(Z) > 0) {
      meanX <- meanX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta_params[3:(2 + length(Z))], ncol = 1))
    }
    scaleX <- meanX / shapeX
    if (any(c(shapeX, scaleX) <= 0)) { return(99999)}
    # ------------------- Shape and scale of Gamma > 0
  } else if (distX == "inverse-gaussian") {
    # Shape and mean of inverse-Gaussian both > 0 ----
    shapeX <- eta_params[1]
    meanX <- eta_params[2]
    if (length(Z) > 0) {
      meanX <- meanX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta_params[3:(2 + length(Z))], ncol = 1))
    }
    if (any(c(shapeX, meanX) <= 0)) { return(99999)}
    # --------------------------------- Get parameters
  } else if (distX == "weibull") {
    # Shape and scale of Weibull both > 0 ------------
    shapeX <- eta_params[1]
    scaleX <- eta_params[2]
    if (length(Z) > 0) {
      eta1 <- eta_params[3:(2 + length(Z))]
      scaleX <- scaleX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
    }
    if (any(c(shapeX, scaleX) <= 0)) { return(99999)}
  } else if (distX %in% c("exponential", "poisson")) {
    # Rate of Exponential or Poisson > 0 -------------
    rateX <- eta_params[1]
    if (length(Z) > 0) {
      eta1 <- eta_params[2:(1 + length(Z))]
      if (length(eta1) == 1) {
        rateX <- rateX + eta1 * uncens_data[, Z]
      } else {
        rateX <- rateX + as.numeric(data.matrix(uncens_data[, Z]) %*% matrix(data = eta1, ncol = 1))
      }
    }
    if (any(rateX <= 0)) { return (99999) }
  }
  pXgivZ <- calc_pXgivZ(x = uncens_data[, X], z = uncens_data[, Z], distX = distX, eta_params = eta_params)

  ####################################################
  # Calculate joint density P(Y,X,Z) #################
  ####################################################
  uncens_data <- cbind(uncens_data, jointP = 1)
  uncens_data[, "jointP"] <- pYgivXZ * pXgivZ
  #uncens_data <- data.frame(cbind(uncens_data, jointP = pYgivXZ * pXgivZ))

  ####################################################
  # Calculate the log-likelihood #####################
  ####################################################
  # Log-likelihood contribution of uncensored X ------
  ll <- sum(log(uncens_data[, "jointP"]))
  # ------ Log-likelihood contribution of uncensored X

  if (nrow(cens_data) > 0) {
    # Log-likelihood contribution of censored X --------
    joint_dens <- function(x, Yi, Zi) {
      ####################################################
      # Analysis model P(Y|X,Z) ##########################
      ####################################################
      pYgivXZ <- calc_pYgivXandZ(y = Yi, x = x, z = Zi, distY = distY, beta_params = beta_params)

      ####################################################
      # Predictor model P(X|Z) ###########################
      ####################################################
      pXgivZ <- calc_pXgivZ(x = x, z = Zi, distX = distX, eta_params = eta_params)

      ####################################################
      # Joint density P(Y,X,Z) ###########################
      ####################################################
      return(pYgivXZ * pXgivZ)
    }
    integrate_joint_dens <- function(data_row) {
      data_row <- data.frame(t(data_row))
      return(
        tryCatch(expr = integrate(f = joint_dens, lower = data_row[, W], upper = Inf, subdivisions = subdivisions,
                                  Yi = data_row[, Y], Zi = data_row[, Z])$value,
                 error = function(err) {0})
      )
    }
    integral <- apply(X = cens_data, MARGIN = 1, FUN = integrate_joint_dens)
    log_integral <- log(integral)
    log_integral[log_integral == -Inf] <- 0
    ll <- ll + sum(log_integral)
    # -------- Log-likelihood contribution of censored X
  }

  # Return (-1) x log-likelihood for use with nlm() --
  return(- ll)
  # -- Return (-1) x log-likelihood for use with nlm()
}
