#' @export
cc_loglik = function(params, dataObj, returnSum = TRUE) {
  ####################################################
  # Create log-likelihood object #####################
  ####################################################
  ## Add all parameters to the data object
  loglikObj = dataObj
  loglikObj$params = params

  ## Add separate beta and eta parameters
  loglikObj$beta_params = get_beta_params(loglikObj)
  loglikObj$eta_params = params[-c(1:length(loglikObj$beta_params))]

  ####################################################
  # Joint conditional density P(Y,X|Z) ###############
  ####################################################
  pYXgivZ = calc_pYXgivZ(x = with(loglikObj, data[, W]),
                         y = with(loglikObj, data[, Y]),
                         z = with(loglikObj, data[, Z]),
                         object = loglikObj)

  ####################################################
  # Likelihood #######################################
  ####################################################
  if (any(is.na(pYXgivZ))) {
    # If params are out of domain, calc_pYXandZ returns NA
    ## And the log-likelihood needs to be arbitrarily "huge"
    return(1E8)
  } else {
    # Replace pYXandZ = 0 with 1 so that log(pYXandZ) = 0
    ## Otherwise log(pYXandZ) = -Inf
    pYXgivZ[pYXgivZ == 0] = 1

    if (returnSum) {
      # Sum over all observations
      ll = sum(log(pYXgivZ))

      # Return (-1) x log-likelihood for use with nlm() --
      return(- ll)
    } else {
      ll = log(pYXgivZ)

      # Return individual log-likelihood contributions
      return(ll)
    }
  }
}

#' @export
cc_loglik_outcome = function(params, dataObj, returnSum = TRUE) {
  ####################################################
  # Create log-likelihood object #####################
  ####################################################
  ## Add all parameters to the data object
  loglikObj = dataObj
  loglikObj$params = params

  ## Add separate beta and eta parameters
  loglikObj$beta_params = params

  ####################################################
  # Conditional density P(Y|X,Z) #####################
  ####################################################
  pYgivXZ = calc_pYgivXZ(x = with(loglikObj, data[, W]),
                         y = with(loglikObj, data[, Y]),
                         z = with(loglikObj, data[, Z]),
                         object = loglikObj)

  ####################################################
  # Likelihood #######################################
  ####################################################
  if (any(is.na(pYgivXZ))) {
    # If params are out of domain, calc_pYXandZ returns NA
    ## And the log-likelihood needs to be arbitrarily "huge"
    return(1E8)
  } else {
    # Replace pYXandZ = 0 with 1 so that log(pYXandZ) = 0
    ## Otherwise log(pYXandZ) = -Inf
    pYgivXZ[pYgivXZ == 0] = 1

    if (returnSum) {
      # Sum over all observations
      ll = sum(log(pYgivXZ))

      # Return (-1) x log-likelihood for use with nlm() --
      return(- ll)
    } else {
      ll = log(pYgivXZ)

      # Return individual log-likelihood contributions
      return(ll)
    }
  }
}

cc_loglik_covariate = function(params, dataObj, returnSum = TRUE) {
  ####################################################
  # Create log-likelihood object #####################
  ####################################################
  ## Add all parameters to the data object
  loglikObj = dataObj
  loglikObj$params = params

  ## Add separate beta and eta parameters
  loglikObj$eta_params = params

  ####################################################
  # Conditional density P(X|Z) #######################
  ####################################################
  pXgivZ = calc_pXgivZ(x = with(loglikObj, data[, W]),
                       z = with(loglikObj, data[, Z]),
                       object = loglikObj)

  ####################################################
  # Likelihood #######################################
  ####################################################
  if (any(is.na(pXgivZ))) {
    # If params are out of domain, calc_pXgivZ returns NA
    ## And the log-likelihood needs to be arbitrarily "huge"
    return(1E8)
  } else {
    # Replace pYXandZ = 0 with 1 so that log(pYXandZ) = 0
    ## Otherwise log(pYXandZ) = -Inf
    pXgivZ[pXgivZ == 0] = 1

    if (returnSum) {
      # Sum over all observations
      ll = sum(log(pXgivZ))

      # Return (-1) x log-likelihood for use with nlm() --
      return(- ll)
    } else {
      ll = log(pXgivZ)

      # Return individual log-likelihood contributions
      return(ll)
    }
  }
}
