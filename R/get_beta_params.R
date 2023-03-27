# Get parameters for model of Y|X,Z
get_beta_params = function(object) {
  UseMethod("get_beta_params")
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
