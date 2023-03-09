calc_deriv2_loglik = function(params = NULL, dataObj, subdivisions) {
  # If no params supplied, take them from dataObj
  if (is.null(params)) {
    params = dataObj$params
  }

  # Save constants
  p = length(params)
  eps = params * (10 ^ (- 4))

  # Create matrix to save derivatives in
  d_theta = matrix(data = 0,
                   nrow = nrow(data),
                   ncol = p ^ 2)
  start_col = 1

  # Calculate derivatives with respect to each pair of parameters
  for (j in 1:p) {
    # Create jth euclidean vector
    ej = matrix(data = 0, nrow = p, ncol = 1)
    ej[j] = 1

    # Perturb 1: params0 = params - eps
    params0 = matrix(data = params - eps * ej,
                     nrow = p,
                     ncol = 1)
    d_l0 = calc_deriv_loglik(params = params0,
                             dataObj = dataObj,
                             subdivisions = subdivisions)
      #calc_deriv_loglik(params = params0, Y = Y, X = X, W = W, D = D, Z = Z, subdivisions = subdivisions, distY = distY, distX = distX, data = data)

    # Perturb 2: params1 = params + eps
    params1 = matrix(data = params + eps * ej,
                     nrow = p,
                     ncol = 1)
    d_l1 = calc_deriv_loglik(params = params1,
                             dataObj = dataObj,
                             subdivisions = subdivisions)
      #calc_deriv_loglik(params = params1, Y = Y, X = X, W = W, D = D, Z = Z, subdivisions = subdivisions, distY = distY, distX = distX, data = data)

    # Estimate deriv = (l'(beta + eps) - l'(beta - eps)) / (2 * eps)
    d_theta[, start_col:(start_col + (p - 1))] = (d_l1 - d_l0) / (2 * eps[j])
    start_col = start_col + p
  }

  return(d_theta)
}
