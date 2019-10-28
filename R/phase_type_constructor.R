#' Constructor function for the class \code{phase_type}.
#'
#' @param type \code{'T_MRCA'}, \code{'T_Total'} or \code{NULL} (default).
#' @param n integer larger than 1 or \code{NULL} (default).
#' @param subint_mat matrix or \code{NULL} (default).
#' @param init_probs vector, a one-row matrix or \code{NULL} (default).
#'
#' @export


phase_type <- function(type = NULL, n = NULL, subint_mat = NULL, init_probs = NULL) {
  if (is.null(type)) {
    if (is.null(subint_mat)) {
      stop('Unable to construct the phase-type distribution. Please provide either the type or the subintensity matrix.')
    }
    else if (is.matrix(subint_mat)) {
      if (is.null(init_probs)) {
        init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
        warning('The initial probability vector is automatically generated.')
      }
      else if ((is.vector(init_probs) & is.atomic(init_probs)) | is.matrix(init_probs)) {
        if (nrow(subint_mat) == length(init_probs)) {
          init_probs <- matrix(init_probs, nrow = 1)
        } else {
          stop('The length of the initial probabilities does not match the size of the subintensity matrix.')
        }
      } else {
        stop('The initial probabilities must be a a matrix with one row or a vector.')
      }
    } else {
      stop('The subintensity matrix must be a matrix.')
    }
  } else if (type == "T_MRCA" & n%%1==0 & n>1) {
    subint_mat = matrix(c(0), nrow = n-1, ncol = n-1)
    for (i in 1:n-1) {
      subint_mat[i, i] = - (n - i+1)*(n - i)/2
      if (i < n-1) {
        subint_mat[i, i+1] = -subint_mat[i, i]
      }
    }
    init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
  } else if (type == "T_Total" & n%%1==0 & n>1) {
    subint_mat = matrix(c(0), nrow = n-1, ncol = n-1)
    for (i in 1:n-1) {
      subint_mat[i, i] = - 0.5 * (n - i)
      if (i < n-1) {
        subint_mat[i, i+1] = -subint_mat[i, i]
      }
    }
    init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
  } else {
    stop('Please provide a valid type')
  }
  value <- list(subint_mat = subint_mat, init_probs = init_probs)
  attr(value, "class") <- "phase_type"
  value
}

#' @export

summary.phase_type <- function(obj) {
  cat('Subintensity matrix:\n')
  print(obj$subint_mat)
  cat('\nInitial probabilities:\n')
  print(obj$init_probs)
}
