#' Create complete dataset
#'
#' @param Y Name of column variable.
#' @param X Name of censored predictor variable.
#' @param W Name of observed (i.e., censored) version of \code{X}.
#' @param D Name of event indicator, defined to be = 1 if \code{X} was uncensored.
#' @param Z (Optional) name(s) of additional fully observed covariates. Default is \code{NULL}.
#' @param partX Size of partition of unobserved \code{X} for censored subjects. Default is \code{50}.
#' @param data A dataframe containing at least columns \code{Y}, \code{X}, \code{W}, \code{D}, and \code{Z}.
#'
#' @return A dataframe containing the same columns as \code{data} but with duplicate rows for each censored subject with different potential \code{X} values.
#'
completeCensRd <- function(Y, X, W, D, Z = NULL, partX = 50, data) {
  # Reorder with uncensored people first ------------------------
  data <- data[order(data[, D], decreasing = TRUE), ]
  n <- nrow(data) # < number of (total) subjects > --------------
  data <- cbind(id = 1:n, data) # Create subject IDs ------------
  n1 <- sum(data[, D]) # < number of uncensored subjects > ------
  # Take uncensored observations "as is" ------------------------
  complete_data_v <- data[1:n1, ]
  ## Pair censored observations with each possible xk based on --
  ### their partition from Wi, ..., max(Xi) ---------------------
  try_X <- unlist(sapply(X = data[(n1 + 1):n, W], FUN = function(x) seq(from = x, to = max(data[1:n1, W]), length.out = partX), simplify = F))
  complete_data_uv <- data[rep((n1 + 1):n, each = partX), ]
  complete_data_uv[, X] <- try_X
  complete_data <- rbind(complete_data_v, complete_data_uv)
  ## Exclude censored observations with Wi > max(Xi) ------------
  complete_data <- complete_data[complete_data[, W] < max(complete_data[1:n1, W]), ]
  # Return artificially complete dataset
  return(complete_data)
}
