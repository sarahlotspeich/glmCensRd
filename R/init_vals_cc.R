init_vals_cc = function(dataObj, D, verbose) {
  # Initial parameter values
  ## Naive initial params (for complete-case)
  if (verbose) {
    print("Start with naive initial values:")
  }
  #params0 = init_vals_naive(object = dataObj)
  beta0 = init_valsY(object = dataObj)
  eta0 = init_valsX(object = dataObj)
  params0 = c(beta0, eta0)
  if (verbose) {
    print(params0)
  }

  ## Get complete-case MLEs from naive initial params
  if (verbose) {
    print("Fit complete-case initial values:")
  }
  cc_dataObj = dataObj
  cc_dataObj$data = cc_dataObj$data[cc_dataObj$data[, D] == 1, ]

  cc_mod_outcome = suppressWarnings(
    nlm(f = cc_loglik_outcome,
        p = beta0,
        dataObj = cc_dataObj,
        steptol = steptol,
        iterlim = iterlim,
        hessian = FALSE
    )
  )

  cc_mod_covariate = suppressWarnings(
    nlm(f = cc_loglik_covariate,
        p = eta0,
        dataObj = cc_dataObj,
        steptol = steptol,
        iterlim = iterlim,
        hessian = FALSE
    )
  )

  if (cc_mod_outcome$code <= 2 & cc_mod_covariate$code <= 2 &
      cc_mod_outcome$iterations > 1 & cc_mod_covariate$iterations > 1) {
    params0 = c(cc_mod_outcome$estimate, cc_mod_covariate$estimate)
  }
  if (verbose) {
    print(params0)
  }

  return(params0)
}
