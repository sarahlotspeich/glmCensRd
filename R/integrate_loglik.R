integrate_loglik = function(object, w, y, z, subdivisions) {
  UseMethod("integrate_loglik")
}

integrate_loglik.rightCensRd = function(object, w, y, z, subdivisions) {
  tryCatch(expr =
             integrate(f = calc_pYXgivZ,
                       lower = w,
                       upper = Inf,
                       subdivisions = subdivisions,
                       y = y,
                       z = z,
                       object = object)$value,
           error = function(err) {0}
  )
}

integrate_loglik.leftCensRd = function(object, w, y, z, subdivisions) {
  tryCatch(expr =
             integrate(f = calc_pYXgivZ,
                       lower = -Inf,
                       upper = w,
                       subdivisions = subdivisions,
                       y = y,
                       z = z,
                       object = object)$value,
           error = function(err) {0}
           )
}