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
  pXgivZ = dpois(x = x, lambda = rateX)
  # -------------------------------------- Calculate
  # Check: rate of Poisson > 0 ---------------------
  pXgivZ[rateX <= 0] = NA
  # --------------------- Check: rate of Poisson > 0
  pXgivZ
}
