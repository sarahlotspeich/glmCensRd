# Calculate probabilities/densities from model of Y|X,Z
calc_pYgivXZ = function(object, y, x, z) {
  UseMethod("calc_pYgivXZ")
}

calc_pYgivXZ.normalY = function(object, y, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)
  y = as.numeric(y)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct mean --------------------------------
  meanY = object$beta_params[1] + object$beta_params[2] * x
  if (!is.null(z)) {
    beta2 = matrix(data = object$beta_params[-c(1:2, length(object$beta_params))],
                    ncol = 1)
    meanY = meanY + as.numeric(z %*% beta2)
  }
  ## Estimate sqrt(variance) directly --------------
  sigY = sqrt(object$beta_params[length(object$beta_params)] ^ 2)
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  dnorm(x = y, mean = meanY, sd = sigY)
  # -------------------------------------- Calculate
}

calc_pYgivXZ.lognormalY = function(object, y, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)
  y = as.numeric(y)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct mean --------------------------------
  meanY = object$beta_params[1] + object$beta_params[2] * x
  if (!is.null(z)) {
    beta2 = matrix(data = object$beta_params[-c(1:2, length(object$beta_params))],
                   ncol = 1)
    meanY = meanY + as.numeric(z %*% beta2)
  }
  ## Estimate sqrt(variance) directly --------------
  sigY = sqrt(object$beta_params[length(object$beta_params)] ^ 2)
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  dlnorm(x = y, meanlog = meanY, sdlog = sigY)
  # -------------------------------------- Calculate
}

calc_pYgivXZ.bernoulliY = function(object, y, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)
  y = as.numeric(y)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct mean --------------------------------
  meanY = object$beta_params[1] + object$beta_params[2] * x
  if (!is.null(z)) {
    beta2 = matrix(data = object$beta_params[-c(1:2)],
                   ncol = 1)
    meanY = meanY + as.numeric(z %*% beta2)
  }
  # --------------------------------- Get parameters
  # Calculate --------------------------------------
  pYgivXZ = exp(- (1 - y) * meanY) / (1 + exp(meanY))
  pYgivXZ[y == 0] = 1 - pYgivXZ[y == 0]
  pYgivXZ
  # -------------------------------------- Calculate
}

calc_pYgivXZ.gammaY = function(object, y, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)
  y = as.numeric(y)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Estimate shape directly -----------------------
  shapeY = object$beta_params[1]
  ## Construct mean --------------------------------
  meanY = object$beta_params[2] + object$beta_params[3] * x
  if (!is.null(z)) {
    beta2 = matrix(data = object$beta_params[-c(1:3)],
                    ncol = 1)
    meanY = meanY + as.numeric(z %*% beta2)
  }
  ## Construct scale -------------------------------
  scaleY = meanY / shapeY
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  if (shapeY <= 0) {
    rep(NA, length(y))
  } else {
    pYgivXZ = dgamma(x = y, shape = shapeY, scale = scaleY)
    pYgivXZ[scaleY <= 0] = NA
    pYgivXZ
  }
  # -------------------------------------- Calculate
}

calc_pYgivXZ.weibullY = function(object, y, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)
  y = as.numeric(y)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Estimate shape directly -----------------------
  shapeY = object$beta_params[1]

  ## Construct scale -------------------------------
  scaleY = object$beta_params[2] + object$beta_params[3] * x
  if (!is.null(z)) {
    beta2 = matrix(data = beta_params[-c(1:3)],
                   ncol = 1)
    scaleY = scaleY + as.numeric(z %*% beta2)
  }
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  if (shapeY <= 0) {
    rep(NA, length(y))
  } else {
    pYgivXZ = dweibull(x = y, shape = shapeY, scale = scaleY)
    pYgivXZ[scaleY <= 0] = NA
    pYgivXZ
  }
  # -------------------------------------- Calculate
}

calc_pYgivXZ.exponentialY = function(object, y, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)
  y = as.numeric(y)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct rate  -------------------------------
  rateY = object$beta_params[1] + object$beta_params[2] * x
  if (!is.null(z)) {
    beta2 = matrix(data = beta_params[-c(1:2)],
                   ncol = 1)
    rateY = rateY + as.numeric(z %*% beta2)
  }
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  pYgivXZ = dexp(x = y, rate = rateY)
  # -------------------------------------- Calculate
  # Check: rate of exponential > 0 -----------------
  pYgivXZ[rateY <= 0] = NA
  # ----------------- Check: rate of exponential > 0
  pYgivXZ
}

calc_pYgivXZ.poissonY = function(object, y, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)
  y = as.numeric(y)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct rate  -------------------------------
  rateY = object$beta_params[1] + object$beta_params[2] * x
  if (!is.null(z)) {
    beta2 = matrix(data = beta_params[-c(1:2)],
                    ncol = 1)
    rateY = rateY + as.numeric(z %*% beta2)
  }
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  pYgivXZ = rateY ^ y * exp(- rateY) / factorial(x = y)
  # -------------------------------------- Calculate
  # Check: rate of Poisson > 0 ---------------------
  pYgivXZ[rateY <= 0] = NA
  # --------------------- Check: rate of Poisson > 0
  pYgivXZ
}
