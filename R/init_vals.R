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
  rep(1E-4, 3 + length(object$Z))
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
  c(rep(0, 1 + length(object$Z)), sd(with(object, data[, X]), na.rm = TRUE))
}

init_valsX.lognormalX = function(object) {
  c(rep(0, 1 + length(object$Z)), sd(with(object, data[, X]), na.rm = TRUE))
}

init_valsX.gammaX = function(object) {
  rep(1E-4, 2 + length(object$Z))
}

init_valsX.weibullX = function(object) {
  rep(1E-4, 2 + length(object$Z))
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
init_vals = function(object) {
  c(init_valsY(object), init_valsX(object))
}
