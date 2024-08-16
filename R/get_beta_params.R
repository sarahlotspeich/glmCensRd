# Get parameters for model of Y|X,Z
get_beta_params = function(object) {
  #UseMethod("get_beta_params")
  if (object$distY == "normal") {
    get_beta_params.normalY(object = object)
  } else if (object$distY == "bernoulli") {
    get_beta_params.bernoulliY(object = object)
  } else if (object$distY == "lognormal") {
    get_beta_params.lognormalY(object = object)
  } else if (object$distY == "gamma") {
    get_beta_params.gammaY(object = object)
  } else if (object$distY == "weibull") {
    get_beta_params.weibullY(object = object)
  } else if (object$distY == "exponential") {
    get_beta_params.exponentialY(object = object)
  } else if (object$distY == "poisson") {
    get_beta_params.poissonY(object = object)
  }
}

get_beta_params.normalY = function(object) {
  lengthZ = length(object$Z)
  object$params[1:(3 + lengthZ)]
}

get_beta_params.bernoulliY = function(object) {
  lengthZ = length(object$Z)
  object$params[1:(2 + lengthZ)]
}

get_beta_params.lognormalY = function(object) {
  lengthZ = length(object$Z)
  object$params[1:(3 + lengthZ)]
}

get_beta_params.gammaY = function(object) {
  lengthZ = length(object$Z)
  object$params[1:(3 + lengthZ)]
}

get_beta_params.weibullY = function(object) {
  lengthZ = length(object$Z)
  object$params[1:(3 + lengthZ)]
}

get_beta_params.exponentialY = function(object) {
  lengthZ = length(object$Z)
  object$params[1:(2 + lengthZ)]
}

get_beta_params.poissonY = function(object) {
  lengthZ = length(object$Z)
  object$params[1:(2 + lengthZ)]
}
