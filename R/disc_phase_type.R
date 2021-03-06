#' \code{disc_phase_type} class
#'
#' Description of the class \code{disc_phase_type}, which represents discrete phase-type
#' distributions.
#'
#' \code{disc_phase_type} is the generator function for the discrete phase-type distribution class
#' of the same name, which inherits from \code{list}. The class is generated by supplying a sub-intensity
#' matrix and an optional initial probability vector. If the initial probabilities are not specified,
#' then a they are automatically generated as the first state having a probability of 1 and the rest
#' a probability of 0.
#'
#'
#' @param init_probs vector, a one-row matrix or \code{NULL} (default).
#' @param subint_mat matrix or \code{NULL} (default).
#'
#' @usage disc_phase_type(subint_mat = NULL, init_probs = NULL)
#'
#' @examples
#' subintensity_matrix = matrix(c(0.4, 0, 0, 0.24, 0.4, 0, 0.12, 0.2, 0.5), ncol = 3)
#' disc_phase_type(subintensity_matrix)
#'
#' #---
#'
#' subintensity_matrix = matrix(c(0.4, 0, 0, 0.24, 0.4, 0, 0.12, 0.2, 0.5), ncol = 3)
#' initial_probabilities = c(0.9, 0.1, 0)
#' disc_phase_type(subintensity_matrix, initial_probabilities)
#'
#' @export

disc_phase_type = function(subint_mat = NULL, init_probs = NULL){

  if (is.null(subint_mat)) {
    stop('Unable to construct the discrete phase-type distribution. Please provide either n or the subintensity matrix.')
  }
  else if (is.matrix(subint_mat)) {
    if(!is.numeric(subint_mat) | nrow(subint_mat) != ncol(subint_mat)){
      stop('Subintensity matrix should be a square numerical matrix')
    }
    # rowsums in Subintensity matrix have to be non positive
    else if(sum(subint_mat < 0) != 0){
      stop('The rowsums in subintensity matrix have to be non-positive')
    }

    if (is.null(init_probs)) {
      init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
      warning('The initial probability vector is automatically generated.')
    }
    else if ((is.vector(init_probs) & is.atomic(init_probs)) | is.matrix(init_probs)) {
      if (nrow(subint_mat) == length(init_probs)) {
        init_probs <- matrix(init_probs, nrow = 1)
      }
      else {
        stop('The length of the initial probabilities does not match the size of the subintensity matrix.')
      }
    }
    else {
      stop('The initial probabilities must be a a matrix with one row or a vector.')
    }
  }
  else {
    stop('The subintensity matrix must be a matrix.')
  }
  value = list(subint_mat = subint_mat, init_probs = init_probs, defect = 1-sum(init_probs))
  attr(value, "class") <- "disc_phase_type"
  return(value)
}


#' @export

mean.disc_phase_type <- function(x, ...) {
  mean <- sum(x$init_probs%*%solve(diag(nrow = nrow(x$subint_mat))-x$subint_mat))
  as.numeric(mean+x$defect)
}


#' @export

var.disc_phase_type <- function(obj) {
  variance <- sum(2*obj$init_probs%*%obj$subint_mat%*%solve((diag(nrow = nrow(obj$subint_mat))-obj$subint_mat)%^%2)) +
    mean(obj) -
    mean(obj)^2
  as.numeric(variance)
}


#' @export

summary.disc_phase_type <- function(object, ...) {
  cat('\nSubintensity matrix:\n')
  print(object$subint_mat)
  cat('\nInitial probabilities:\n')
  print(object$init_probs)
  cat('\nDefect:\n')
  print(object$defect)
  cat('\nMean: ', mean(object), '\n', sep = '')
  cat('\nVariance: ', var(object), '\n\n', sep = '')
}
