cc_loglik = function(params, dataObj) {
  ####################################################
  # Create log-likelihood object #####################
  ####################################################
  ## Add all parameters to the data object
  loglikObj = dataObj
  loglikObj$params = params

  ## Add separate beta and eta parameters
  loglikObj$beta_params = get_beta_params(loglikObj)
  loglikObj$eta_params = params[-c(1:length(loglikObj$beta_params))]

  print(loglikObj)

  ####################################################
  # Joint density P(Y,X,Z) ###########################
  ####################################################
  pYXgivZ = calc_pYXgivZ(object = loglikObj,
                         y = loglikObj$data[, loglikObj$Y],
                         x = loglikObj$data[, loglikObj$W],
                         z = loglikObj$data[, loglikObj$Z])

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
    ll = sum(log(pYXgivZ))

    # Return (-1) x log-likelihood for use with nlm() --
    return(- ll)
  }
}