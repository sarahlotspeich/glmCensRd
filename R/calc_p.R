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
  pYgivXZ = exp(- (1 - y) * meanY) / (1 + exp(-meanY))
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

# Calculate probabilities/densities from model of X|Z
calc_pXgivZ = function(object, x, z) {
  UseMethod("calc_pXgivZ")
}

calc_pXgivZ.normalX = function(object, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct mean --------------------------------
  meanX = object$eta_params[1]
  if (!is.null(z)) {
    eta1 = matrix(data = object$eta_params[-c(1, length(object$eta_params))],
                  ncol = 1)
    meanX = meanX + as.numeric(z %*% eta1)
  }
  ## Estimate sqrt(variance) directly --------------
  sigX = sqrt(object$eta_params[length(object$eta_params)] ^ 2)
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  dnorm(x = x, mean = meanX, sd = sigX)
  # -------------------------------------- Calculate
}

calc_pXgivZ.lognormalX = function(object, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct mean --------------------------------
  meanX = object$eta_params[1]
  if (!is.null(z)) {
    eta1 = matrix(data = object$eta_params[-c(1, length(object$eta_params))],
                  ncol = 1)
    meanX = meanX + as.numeric(z %*% eta1)
  }
  ## Estimate sqrt(variance) directly --------------
  sigX = sqrt(object$eta_params[length(object$eta_params)] ^ 2)
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  dlnorm(x = x, meanlog = meanX, sdlog = sigX)
  # -------------------------------------- Calculate
}

calc_pXgivZ.gammaX = function(object, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Estimate shape directly -----------------------
  shapeX = object$eta_params[1]
  ## Construct mean --------------------------------
  meanX = object$eta_params[2]
  if (!is.null(z)) {
    eta1 = matrix(data = object$eta_params[-c(1:2)],
                  ncol = 1)
    meanX = meanX + as.numeric(z %*% eta1)
  }
  ## Construct scale -------------------------------
  scaleX = meanX / shapeX
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  if (shapeX <= 0) {
    rep(NA, length(x))
  } else {
    pXgivZ = dgamma(x = x, shape = shapeX, scale = scaleX)
    pXgivZ[scaleX <= 0] = NA
    pXgivZ
  }
  # -------------------------------------- Calculate
}

calc_pXgivZ.weibullX = function(object, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Estimate shape directly -----------------------
  shapeX = object$eta_params[1]

  ## Construct scale -------------------------------
  scaleX = object$eta_params[2]
  if (!is.null(z)) {
    eta1 = matrix(data = object$eta_params[-c(1:2)],
                  ncol = 1)
    scaleX = scaleX + as.numeric(z %*% eta1)
  }
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  if (shapeX <= 0) {
    rep(NA, length(x))
  } else {
    pXgivZ = dweibull(x = x, shape = shapeX, scale = scaleX)
    pXgivZ[scaleX <= 0] = NA
    pXgivZ
  }
  # -------------------------------------- Calculate
}

calc_pXgivZ.exponentialX = function(object, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct rate  -------------------------------
  rateX = object$eta_params[1]
  if (!is.null(z)) {
    eta1 = matrix(data = object$eta_params[-c(1)],
                  ncol = 1)
    rateX = rateX + as.numeric(z %*% eta1)
  }
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  pXgivZ = dexp(x = x, rate = rateX)
  # -------------------------------------- Calculate
  # Check: rate of exponential > 0 -----------------
  pXgivZ[rateX <= 0] = NA
  # ----------------- Check: rate of exponential > 0
  pXgivZ
}

calc_pXgivZ.poissonX = function(object, x, z) {
  # Treat y and x as numeric vectors (not dataframe, matrix, etc)
  x = as.numeric(x)

  # Treat z as a matrix
  z = data.matrix(z)

  # Get parameters ---------------------------------
  ## Construct rate  -------------------------------
  rateX = object$eta_params[1]
  if (!is.null(z)) {
    eta1 = matrix(data = object$eta_params[-c(1)],
                  ncol = 1)
    rateX = rateX + as.numeric(z %*% eta1)
  }
  # --------------------------------- Get parameters

  # Calculate --------------------------------------
  pXgivZ = rateX ^ x * exp(- rateX) / factorial(x = x)
  # -------------------------------------- Calculate
  # Check: rate of Poisson > 0 ---------------------
  pXgivZ[rateX <= 0] = NA
  # --------------------- Check: rate of Poisson > 0
  pXgivZ
}

# Calculate probabilities/densities from model of Y,X|Z
calc_pYXgivZ = function(x, y, z, object) {
  calc_pYgivXZ(object = object, y = y, x = x, z = z) *
    calc_pXgivZ(object = object, x = x, z = z)
}
