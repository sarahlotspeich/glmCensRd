# Get initial values for model of Y|X,Z
init_valsY = function(object) {
  UseMethod("init_valsY")
}

init_valsY.normalY = function(object) {
  c(rep(0, 2 + length(object$Z)), sd(with(object, data[, Y])))
}

init_valsY.lognormalY = function(object) {
  c(rep(0, 2 + length(object$Z)), sd(with(object, data[, Y])))
}

init_valsY.bernoulliY = function(object) {
  rep(0, 2 + length(object$Z))
}

init_valsY.gammaY = function(object) {
  rep(1E-4, 3 + length(object$Z))
}

init_valsY.weibullY = function(object) {
  c(1, rep(1E-4, 2 + length(object$Z)))
}

init_valsY.inversegaussianY = function(object) {
  rep(1E-4, 3 + length(object$Z))
}

init_valsY.exponentialY = function(object) {
  c(1E-4, rep(0, 1 + length(object$Z)))
}

init_valsY.poissonY = function(object) {
  c(1E-4, rep(0, 1 + length(object$Z)))
}

# Get initial values for model of X|Z
init_valsX = function(object) {
  UseMethod("init_valsX")
}

init_valsX.normalX = function(object) {
  c(rep(0, 1 + length(object$Z)), sd(with(object, data[, "X"]), na.rm = TRUE))
}

init_valsX.lognormalX = function(object) {
  c(rep(0, 1 + length(object$Z)), sd(with(object, data[, "X"]), na.rm = TRUE))
}

init_valsX.gammaX = function(object) {
  rep(1E-4, 2 + length(object$Z))
}

init_valsX.weibullX = function(object) {
  c(1, rep(1E-4, 1 + length(object$Z)))
}

init_valsX.inversegaussianX = function(object) {
  rep(1E-4, 2 + length(object$Z))
}

init_valsX.exponentialX = function(object) {
  c(1E-4, rep(0, length(object$Z)))
}

init_valsX.poissonX = function(object) {
  c(1E-4, rep(0, length(object$Z)))
}

# Get initial values for both models
init_vals_naive = function(object) {
  c(init_valsY(object), init_valsX(object))
}

init_vals_cc = function(dataObj, D, steptol, iterlim, verbose) {
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
