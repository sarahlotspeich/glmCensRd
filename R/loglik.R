#' @export
loglik = function(params, dataObj, subdivisions){ #params, Y, X, W, D, Z = NULL, data, distY = "normal", distX = "normal", cens = "right", subdivisions = 100) {
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
  ll = - cc_loglik(params = params,
                   dataObj = uncens_dataObj)

  ####################################################
  # Log-likelihood of censored observations ##########
  ####################################################
  #if (nrow(cens_dataObj$data) > 0) {
    integrate_pYXgivZ = function(data_row) {
      Wi = as.numeric(data_row[cens_dataObj$W])
      Yi = data_row[cens_dataObj$Y]
      Zi = data_row[cens_dataObj$Z]
      return(
        #tryCatch(expr =
        integrate(f = calc_pYXgivZ,
                                  lower = ifelse(test = cens_dataObj$rightCens,
                                                 yes = Wi,
                                                 no = -Inf),
                                  upper = ifelse(test = cens_dataObj$rightCens,
                                                 yes = Inf,
                                                 no = Wi),
                                  subdivisions = subdivisions,
                                  y = Yi,
                                  z = Zi,
                                  object = cens_dataObj)$value#,
                 #error = function(err) {0})
      )
    }
    int_pYXgivZ_cens = apply(X = cens_dataObj$data,
                             MARGIN = 1,
                             FUN = integrate_pYXgivZ)
    log_int_pYXgivZ_cens = log(int_pYXgivZ_cens)
    log_int_pYXgivZ_cens[log_int_pYXgivZ_cens == -Inf] = 0
    ll = ll + sum(log_int_pYXgivZ_cens)
  #}

  # Return (-1) x log-likelihood for use with nlm() --
  return(- ll)
  # -- Return (-1) x log-likelihood for use with nlm()
}
