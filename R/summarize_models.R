summarize_models = function(est, se, rob_se, X, Z = NULL, distY = "normal", distX = "normal") {
  ####################################################
  # Analysis model P(Y|X,Z) ##########################
  ####################################################
  if (distY == "normal") {
    # Create coefficients dataframe for outcome model
    ## Mean parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_mean_est = est[1:dim_beta]
    modY_mean_se = se[1:dim_beta]
    modY_mean_rse = rob_se[1:dim_beta]
    modY_mean = data.frame(coeff = modY_mean_est,
                           se = modY_mean_se,
                           robse = modY_mean_rse)
    rownames(modY_mean) = c("(Intercept)", X, Z)
    
    ## Error variance parameter (estimated directly)
    modY_sigma2 = est[dim_beta + 1] ^ 2
    
    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                mean = modY_mean,
                sigma2 = modY_sigma2)
  } else if (distY == "binomial") {
    # Create coefficients dataframe for outcome model
    ## Mean parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_mean_est = est[1:dim_beta]
    modY_mean_se = se[1:dim_beta]
    modY_mean_rse = rob_se[1:dim_beta]
    modY_mean = data.frame(coeff = modY_mean_est,
                           se = modY_mean_se,
                           robse = modY_mean_rse)
    rownames(modY_mean) = c("(Intercept)", X, Z)
    
    # Construct contents of "covariate_model" slot
    modY = list(distY = distY,
                mean = modY_mean)
  } else if (distY %in% c("gamma", "inverse-gaussian")) {
    # Create coefficients dataframe for covariate model
    ## Shape parameter (estimated directly)
    modY_shape_est = est[1]
    modY_shape_se = se[1]
    modY_shape_rse = rob_se[1]
    modY_shape = data.frame(coeff = modY_shape_est,
                            se = modY_shape_se,
                            robse = modY_shape_rse)
    rownames(modY_shape) = c("(Intercept)")
    est = est[-1]
    se = se[-1]
    rob_se = rob_se[-1]
    
    ## Mean parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_mean_est = est[1:dim_beta]
    modY_mean_se = se[1:dim_beta]
    modY_mean_rse = rob_se[1:dim_beta]
    modY_mean = data.frame(coeff = modY_mean_est,
                           se = modY_mean_se,
                           robse = modY_mean_rse)
    rownames(modY_mean) = c("(Intercept)", X, Z)
    
    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                mean = modY_mean,
                shape = modY_shape)
  } else if (distY == "weibull") {
    ## Shape parameter (estimated directly)
    modY_shape_est = est[1]
    modY_shape_se = se[1]
    modY_shape_rse = rob_se[1]
    modY_shape = data.frame(coeff = modY_shape_est,
                            se = modY_shape_se,
                            robse = modY_shape_rse)
    rownames(modY_shape) = c("(Intercept)")
    est = est[-1]
    se = se[-1]
    rob_se = rob_se[-1]
    
    ## Scale parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_scale_est = est[1:dim_beta]
    modY_scale_se = se[1:dim_beta]
    modY_scale_rse = rob_se[1:dim_beta]
    modY_scale = data.frame(coeff = modY_scale_est,
                            se = modY_scale_se,
                            robse = modY_scale_rse)
    rownames(modY_scale) = c("(Intercept)", X, Z)
    
    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                scale = modY_scale,
                shape = modY_shape)
  } else if (distY %in% c("exponential", "poisson")) {
    # Create coefficients dataframe for outcome model
    ## Rate parameter (linear function of Z)
    dim_beta = length(Z) + 2
    modY_rate_est = est[1:dim_beta]
    modY_rate_se = se[1:dim_beta]
    modY_rate_rse = rob_se[1:dim_beta]
    modY_rate = data.frame(coeff = modY_rate_est,
                           se = modY_rate_se,
                           robse = modY_rate_rse)
    rownames(modY_rate) = c("(Intercept)", X, Z)
    
    # Construct contents of "outcome_model" slot
    modY = list(distY = distY,
                rate = modY_rate)
  }
  # Remove outcome model parameters/standard errors
  est = est[-c(1:(dim_beta + 1))]
  se = se[-c(1:(dim_beta + 1))]
  rob_se = rob_se[-c(1:(dim_beta + 1))]
  
  ####################################################
  # covariate model P(X|Z) ###########################
  ####################################################
  if (distX %in% c("normal", "log-normal")) {
    # Create coefficients dataframe for covariate model
    ## Mean parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_mean_est = est[1:dim_eta]
    modX_mean_se = se[1:dim_eta]
    modX_mean_rse = rob_se[1:dim_eta]
    modX_mean = data.frame(coeff = modX_mean_est,
                           se = modX_mean_se,
                           robse = modX_mean_rse)
    rownames(modX_mean) = c("(Intercept)", Z)
    
    ## Error variance parameter (estimated directly)
    modX_sigma2 = est[dim_eta + 1] ^ 2
    
    # Construct contents of "covariate_model" slot
    modX = list(distX = distX,
                mean = modX_mean,
                sigma2 = modX_sigma2)
  } else if (distX %in% c("gamma", "inverse-gaussian")) {
    # Create coefficients dataframe for covariate model
    ## Shape parameter (estimated directly)
    modX_shape_est = est[1]
    modX_shape_se = se[1]
    modX_shape_rse = rob_se[1]
    modX_shape = data.frame(coeff = modX_shape_est,
                            se = modX_shape_se,
                            robse = modX_shape_rse)
    rownames(modX_shape) = c("(Intercept)")
    est = est[-1]
    se = se[-1]
    rob_se = rob_se[-1]
    
    ## Mean parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_mean_est = est[1:dim_eta]
    modX_mean_se = se[1:dim_eta]
    modX_mean_rse = rob_se[1:dim_eta]
    modX_mean = data.frame(coeff = modX_mean_est,
                           se = modX_mean_se,
                           robse = modX_mean_rse)
    rownames(modX_mean) = c("(Intercept)", Z)
    
    # Construct contents of "covariate_model" slot
    modX = list(distX = distX,
                mean = modX_mean,
                shape = modX_shape)
  } else if (distX == "weibull") {
    ## Shape parameter (estimated directly)
    modX_shape_est = est[1]
    modX_shape_se = se[1]
    modX_shape_rse = rob_se[1]
    modX_shape = data.frame(coeff = modX_shape_est,
                            se = modX_shape_se,
                            robse = modX_shape_rse)
    rownames(modX_shape) = c("(Intercept)")
    est = est[-1]
    se = se[-1]
    rob_se = rob_se[-1]
    
    ## Scale parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_scale_est = est[1:dim_eta]
    modX_scale_se = se[1:dim_eta]
    modX_scale_rse = rob_se[1:dim_eta]
    modX_scale = data.frame(coeff = modX_scale_est,
                            se = modX_scale_se,
                            robse = modX_scale_rse)
    rownames(modX_scale) = c("(Intercept)", Z)
    
    # Construct contents of "covariate_model" slot
    modX = list(distX = distX,
                scale = modX_scale,
                shape = modX_shape)
  } else if (distX %in% c("exponential", "poisson")) {
    ## Rate parameter (linear function of Z)
    dim_eta = length(Z) + 1
    modX_rate_est = est[1:dim_eta]
    modX_rate_se = se[1:dim_eta]
    modX_rate_rse = rob_se[1:dim_eta]
    modX_rate = data.frame(coeff = modX_rate_est,
                           se = modX_rate_se,
                           robse = modX_rate_rse)
    rownames(modX_rate) = c("(Intercept)", Z)
    
    # Construct contents of "covariate_model" slot
    modX = list(distX = distX,
                rate = modX_rate)
  }
  
  return(list(modY = modY, 
              modX = modX))
}
