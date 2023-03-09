calc_deriv_loglik = function(params = NULL, dataObj, subdivisions) {
  # If no params supplied, take them from dataObj
  if (is.null(params)) {
    params = dataObj$params
  }

  # Save constants
  p = length(params)
  eps = params * (10 ^ (- 4))

  # Create matrix to save derivatives in
  d_theta = matrix(data = 0,
                   nrow = nrow(dataObj$data),
                   ncol = p)

  # Calculate derivative with respect to each parameter
  for (j in 1:p) {
    # Create jth euclidean vector
    ej = matrix(data = 0,
                nrow = p,
                ncol = 1)
    ej[j] = 1

    # Perturb 1: params0 = params - eps
    params0 = matrix(data = params - eps * ej,
                     nrow = p,
                     ncol = 1)
    l0 = loglik(params = params0,
                dataObj = dataObj,
                returnSum = FALSE,
                subdivisions = subdivisions)

    # l0 = calc_indiv_loglik(params = params0,
    #                        Y = Y,
    #                        X = X,
    #                        W = W,
    #                        D = D,
    #                        Z = Z,
    #                        subdivisions = subdivisions,
    #                        distY = distY,
    #                        distX = distX,
    #                        data = data)

    # Perturb 2: params1 = params + eps
    params1 = matrix(data = params + eps * ej,
                     nrow = p,
                     ncol = 1)
    l1 = loglik(params = params1,
                dataObj = dataObj,
                returnSum = FALSE,
                subdivisions = subdivisions)

    # l1 = calc_indiv_loglik(params = params1,
    #                        Y = Y,
    #                        X = X,
    #                        W = W,
    #                        D = D,
    #                        Z = Z,
    #                        subdivisions = subdivisions,
    #                        distY = distY,
    #                        distX = distX,
    #                        data = data)

    # Estimate deriv = (l(beta + eps) - l(beta - eps)) / (2 * eps)
    d_theta[, j] = (l1 - l0) / (2 * eps[j])
  }

  return(d_theta)
}
