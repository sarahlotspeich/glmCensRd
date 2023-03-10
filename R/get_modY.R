# Calculate probabilities/densities from model of Y|X,Z
get_modY = function(object) {
  UseMethod("get_modY")
}

get_modY.normalY = function(object) {
  # Create coefficients dataframe for outcome model
  ## Mean parameter (linear function of Z)
  dim_beta = length(object$Z) + 2
  mean_est = object$params[1:dim_beta]
  mean_se = object$se[1:dim_beta]
  mean_rse = object$rob_se[1:dim_beta]
  mean = data.frame(coeff = mean_est,
                         se = mean_se,
                         robse = mean_rse)
  rownames(mean) = c("(Intercept)", "X", object$Z)

  ## Error variance parameter (estimated directly)
  sigma2 = object$params[dim_beta + 1] ^ 2

  # Construct contents of "outcome_model" slot
  modY = list(mean = mean,
              sigma2 = sigma2)
  modY
}

get_modY.bernoulliY = function(object) {
  # Create coefficients dataframe for outcome model
  ## Mean parameter (linear function of Z)
  dim_beta = length(object$Z) + 2
  mean_est = object$params[1:dim_beta]
  mean_se = object$se[1:dim_beta]
  mean_rse = object$rob_se[1:dim_beta]
  mean = data.frame(coeff = mean_est,
                    se = mean_se,
                    robse = mean_rse)
  rownames(mean) = c("(Intercept)", "X", object$Z)

  # Construct contents of "outcome_model" slot
  modY = list(mean = mean)
  modY
}

get_modY.lognormalY = function(object) {
  # Create coefficients dataframe for outcome model
  ## Mean parameter (linear function of Z)
  dim_beta = length(object$Z) + 2
  mean_est = object$params[1:dim_beta]
  mean_se = object$se[1:dim_beta]
  mean_rse = object$rob_se[1:dim_beta]
  mean = data.frame(coeff = mean_est,
                    se = mean_se,
                    robse = mean_rse)
  rownames(mean) = c("(Intercept)", "X", object$Z)

  ## Error variance parameter (estimated directly)
  sigma2 = object$params[dim_beta + 1] ^ 2

  # Construct contents of "outcome_model" slot
  modY = list(mean = mean,
              sigma2 = sigma2)
  modY
}

get_modY.gammaY = function(object) {
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
  dim_beta = length(object$Z) + 2
  mean_est = object$params[1:dim_beta]
  mean_se = object$se[1:dim_beta]
  mean_rse = object$rob_se[1:dim_beta]
  mean = data.frame(coeff = mean_est,
                    se = mean_se,
                    robse = mean_rse)
  rownames(mean) = c("(Intercept)", "X", object$Z)

  # Construct contents of "outcome_model" slot
  modY = list(mean = mean,
              shape = shape)
  modY
}

get_modY.weibullY = function(object) {
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
  dim_beta = length(object$Z) + 2
  scale_est = param_est[1:dim_beta]
  scale_se = param_se[1:dim_beta]
  scale_rse = param_rob_se[1:dim_beta]
  scale = data.frame(coeff = scale_est,
                          se = scale_se,
                          robse = scale_rse)
  rownames(scale) = c("(Intercept)", "X", object$Z)

  # Construct contents of "outcome_model" slot
  modY = list(scale = scale,
              shape = shape)
  modY
}

get_modY.exponentialY = function(object) {
  # Create coefficients dataframe for outcome model
  ## Rate parameter (linear function of Z)
  rate_est = object$params
  rate_se = object$se
  rate_rse = object$rob_se
  rate = data.frame(coeff = rate_est,
                    se = rate_se,
                    robse = rate_rse)
  rownames(rate) = c("(Intercept)", "X", object$Z)

  # Construct contents of "outcome_model" slot
  modY = list(rate = rate)
  modY
}

get_modY.poissonY = function(object) {
  # Create coefficients dataframe for outcome model
  ## Rate parameter (linear function of Z)
  rate_est = object$params
  rate_se = object$se
  rate_rse = object$rob_se
  rate = data.frame(coeff = rate_est,
                    se = rate_se,
                    robse = rate_rse)
  rownames(rate) = c("(Intercept)", "X", object$Z)

  # Construct contents of "outcome_model" slot
  modY = list(rate = rate)
  modY
}
