conv = function(object) {
  UseMethod("conv")
}

conv.glmCensRd = function(object) {
  object$code
}