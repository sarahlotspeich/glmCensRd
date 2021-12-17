#' Approximate the integral over censored observations in the log-likelihood using quadrature (the trapezoidal rule)
#'
#' @param X Name of censored predictor variable.
#' @param complete_data_cens A dataframe containing artificially complete data for the censored subjects (returned from \code{completeCensRd()})
#'
#' @return A vector containing the approximate integral over the joint density of censored subjects from \code{complete_data_cens}.
#'
integrateCensRd <- function(X, complete_data_cens) {
  # Find x_{1} and x_{m} for each id ---------------
  which_x1 <- unlist(lapply(X = unique(complete_data_cens[, "id"]),
                     FUN = function(x, all_x = complete_data_cens[, "id"]) return(min(which(all_x == x)))))
  which_xm <- unlist(lapply(X = unique(complete_data_cens[, "id"]),
                            FUN = function(x, all_x = complete_data_cens[, "id"]) return(max(which(all_x == x)))))
  # Calculate (x_{k+1} - x_{k}) --------------------
  subtract_X <- complete_data_cens[- which_x1, X] - complete_data_cens[- which_xm, X]
  # Calculate P(Yi, x_{k+1}, Zi) + -----------------
  ## P(Yi, x_{k}, Zi) ------------------------------
  sum_joint <- complete_data_cens[- which_x1, "jointP"] + complete_data_cens[- which_xm, "jointP"]
  # ------------------------------- Create indicator
  # Sum over trapezoids (within id) ----------------
  integral <- 1/2 * rowsum(x = subtract_X * sum_joint, group = complete_data_cens[- which_x1, "id"])
  return(integral)
}
