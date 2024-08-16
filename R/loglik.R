#' @export
loglik = function(params, dataObj, returnSum = TRUE, subdivisions){
  ####################################################
  # Pre-processing ###################################
  ####################################################
  ## Create log-likelihood object
  ### Add all parameters to the data object
  loglikObj = dataObj
  loglikObj$params = params

  ## Add separate beta and eta parameters
  loglikObj$beta_params = get_beta_params(loglikObj)
  loglikObj$eta_params = params[-c(1:length(loglikObj$beta_params))]

  ## Create log-likelihood objects for censored...
  cens_dataObj = loglikObj
  cens_dataObj$data = with(cens_dataObj, data[data[, D] == 0, ])
  ### and uncensored observations
  uncens_dataObj = loglikObj
  uncens_dataObj$data = with(uncens_dataObj, data[data[, D] == 1, ])

  ####################################################
  # Log-likelihood of uncensored observations ########
  ####################################################
  ll = cc_loglik(params = params,
                 dataObj = uncens_dataObj,
                 returnSum = returnSum)
  ## If returnSum, cc_loglik() returns (-1) x log-
  ## likelihood for use with nlm()
  ll = ll * ifelse(test = returnSum,
                   yes = -1, ### cancel this out
                   no = 1)

  ####################################################
  # Log-likelihood of censored observations ##########
  ####################################################
  if (nrow(cens_dataObj$data) > 0) {
    integrate_pYXgivZ = function(data_row) {
      Wi = as.numeric(data_row[cens_dataObj$W])
      Yi = data_row[cens_dataObj$Y]
      Zi = data_row[cens_dataObj$Z]
      return(
        integrate_loglik(object = cens_dataObj,
                         w = Wi,
                         y = Yi,
                         z = Zi,
                         subdivisions = subdivisions)
      )
    }
    int_pYXgivZ_cens = apply(X = cens_dataObj$data,
                             MARGIN = 1,
                             FUN = integrate_pYXgivZ)
    log_int_pYXgivZ_cens = log(int_pYXgivZ_cens)
    log_int_pYXgivZ_cens[log_int_pYXgivZ_cens == -Inf] = 0

    if (returnSum) {
      ll = ll + sum(log_int_pYXgivZ_cens)
    } else {
      ll = c(ll, log_int_pYXgivZ_cens)
    }
  }

  ####################################################
  # Return ###########################################
  ####################################################
  if (returnSum) {
    ## Return (-1) x log-likelihood for use with nlm()
    return(- ll)
  } else {
    ## Return vector of log-likelihoods --------------
    return(ll)
  }
}

integrate_loglik = function(object, w, y, z, subdivisions) {
  #UseMethod("integrate_loglik")
  if (object$rightCens) {
    integrate_loglik.rightCensRd(object = object,
                                 w = w,
                                 y = y,
                                 z = z,
                                 subdivisions = subdivisions)
  } else {
    integrate_loglik.leftCensRd(object = object,
                                w = w,
                                y = y,
                                z = z,
                                subdivisions = subdivisions)
  }
}

integrate_loglik.rightCensRd = function(object, w, y, z, subdivisions) {
  tryCatch(
    expr = integrate(
      f = calc_pYXgivZ,
      lower = w,
      upper = Inf,
      subdivisions = subdivisions,
      y = y,
      z = z,
      object = object)$value,
    error = function(err) {0}
  )
}

integrate_loglik.leftCensRd = function(object, w, y, z, subdivisions) {
  tryCatch(
    expr = integrate(
      f = calc_pYXgivZ,
      lower = -Inf,
      upper = w,
      subdivisions = subdivisions,
      y = y,
      z = z,
      object = object)$value,
    error = function(err) {0}
  )
}
