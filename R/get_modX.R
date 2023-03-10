# Calculate probabilities/densities from model of Y|X,Z
get_modX = function(object) {
  UseMethod("get_modX")
}

get_modX.normalX = function(object) {
  # Create coefficients dataframe for outcome model
  ## Mean parameter (linear function of Z)
  mean_est = object$params[-length(object$params)]
  mean_se = object$se[-length(object$params)]
  mean_rse = object$rob_se[-length(object$params)]
  mean = data.frame(coeff = mean_est,
                    se = mean_se,
                    robse = mean_rse)
  rownames(mean) = c("(Intercept)", object$Z)

  ## Error variance parameter (estimated directly)
  sigma2 = object$params[length(object$params)] ^ 2

  # Construct contents of "covariate_model" slot
  modX = list(mean = mean,
              sigma2 = sigma2)
  class(modX) = class(object)[2]
  modX
}

get_modX.lognormalX = function(object) {
  # Create coefficients dataframe for outcome model
  ## Mean parameter (linear function of Z)
  mean_est = object$params[-length(object$params)]
  mean_se = object$se[-length(object$params)]
  mean_rse = object$rob_se[-length(object$params)]
  mean = data.frame(coeff = mean_est,
                    se = mean_se,
                    robse = mean_rse)
  rownames(mean) = c("(Intercept)", object$Z)

  ## Error variance parameter (estimated directly)
  sigma2 = object$params[length(object$params)] ^ 2

  # Construct contents of "covariate_model" slot
  modX = list(mean = mean,
              sigma2 = sigma2)
  class(modX) = class(object)[2]
  modX
}

get_modX.gammaX = function(object) {
  # Create coefficients dataframe for predictor model
  ## Shape parameter (estimated directly)
  shape_est = object$params[1]
  shape_se = object$se[1]
  shape_rse = object$rob_se[1]
  shape = data.frame(coeff = shape_est,
                     se = shape_se,
                     robse = shape_rse)
  rownames(shape) = c("(Intercept)")

  est = object$params[-1]
  se = object$se[-1]
  rob_se = object$rob_se[-1]

  ## Mean parameter (linear function of Z)
  dim_beta = length(object$Z) + 1
  mean_est = object$params[1:dim_beta]
  mean_se = object$se[1:dim_beta]
  mean_rse = object$rob_se[1:dim_beta]
  mean = data.frame(coeff = mean_est,
                    se = mean_se,
                    robse = mean_rse)
  rownames(mean) = c("(Intercept)", object$Z)

  # Construct contents of "covariate_model" slot
  modX = list(mean = mean,
              shape = shape)
  class(modX) = class(object)[2]
  modX
}

get_modX.weibullX = function(object) {
  ## Shape parameter (estimated directly)
  shape_est = object$params[1]
  shape_se = object$se[1]
  shape_rse = object$rob_se[1]
  shape = data.frame(coeff = shape_est,
                     se = shape_se,
                     robse = shape_rse)
  rownames(shape) = c("(Intercept)")
  param_est = object$params[-1]
  param_se = object$se[-1]
  param_rob_se = object$rob_se[-1]

  ## Scale parameter (linear function of Z)
  dim_beta = length(object$Z) + 1
  scale_est = param_est[1:dim_beta]
  scale_se = param_se[1:dim_beta]
  scale_rse = param_rob_se[1:dim_beta]
  scale = data.frame(coeff = scale_est,
                          se = scale_se,
                          robse = scale_rse)
  rownames(scale) = c("(Intercept)", object$Z)

  # Construct contents of "covariate_model" slot
  modX = list(scale = scale,
              shape = shape)
  class(modX) = class(object)[2]
  modX
}

get_modX.exponentialX = function(object) {
  # Create coefficients dataframe for outcome model
  ## Rate parameter (linear function of Z)
  rate_est = object$params
  rate_se = object$se
  rate_rse = object$rob_se
  rate = data.frame(coeff = rate_est,
                    se = rate_se,
                    robse = rate_rse)
  rownames(rate) = c("(Intercept)", "X", object$Z)

  # Construct contents of "covariate_model" slot
  modX = list(rate = rate)
  class(modX) = class(object)[2]
  modX
}

get_modX.poissonX = function(object) {
  # Create coefficients dataframe for outcome model
  ## Rate parameter (linear function of Z)
  rate_est = object$params
  rate_se = object$se
  rate_rse = object$rob_se
  rate = data.frame(coeff = rate_est,
                    se = rate_se,
                    robse = rate_rse)
  rownames(rate) = c("(Intercept)", "X", object$Z)

  # Construct contents of "covariate_model" slot
  modX = list(rate = rate)
  class(modX) = class(object)[2]
  modX
}
