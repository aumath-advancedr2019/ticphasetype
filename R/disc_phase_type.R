#' \code{disc_phase_type} class
#'
#' Description of the class \code{disc_phase_type}, which represents discrete phase-type
#' distributions.
#'
#' \code{disc_phase_type} is the generator function for the discrete phase-type distribution class
#' of the same name, which inherits from \code{list}. This function can easily produce phase-type
#' representations of common statistics in genomics, such as the total number of segregating sites,
#' the number of i-tons (singletons, doubletons, etc.) and the tail statistic.
#'
#' The discrete phase-type distribution can be generated in two different ways:
#' \itemize{
#'   \item By specifying whether the user wants a  discrete phase-type representation
#'   of the total number of segregating sites (default), the number of i-tons, or the tail
#'   statistic for a certain number of sequences \code{n}.
#'   \item By supplying a user-defined sub-intensity matrix, with optional initial
#'   probabilities.
#' }
#'
#' See examples for further explanation on the usage.
#'
#' @param n integer larger than 1, or \code{NULL} (default).
#' @param itons integer between 1 and n-1, or \code{NULL} (default).
#' @param init_probs vector, a one-row matrix or \code{NULL} (default).
#' @param subint_mat matrix or \code{NULL} (default).
#' @param theta numeric.
#' @param tail_stat logical.
#'
#' @usage disc_phase_type(n = NULL, itons = NULL,
#'             theta = 2, tail_stat = F,
#'             subint_mat = NULL, init_probs = NULL)
#'
#'
#' @examples
#' # Total number of segregating sites
#' disc_phase_type(4)
#' disc_phase_type(5, theta=0.5)
#'
#' # Below examples are not implemented:
#' # Number of singletons
#' disc_ph_example <- disc_phase_type(4, itons=1)
#' mean(disc_ph_example)
#' var(disc_ph_example)
#' summary(disc_ph_example)
#'
#' # Tail statistic
#' disc_ph_tail <- disc_phase_type(4, itons=2, tail_stat=TRUE)
#' summary(disc_ph_tail)
#'
#' @export

disc_phase_type = function(subint_mat = NULL, init_probs = NULL){

  if (is.null(subint_mat)) {
    stop('Unable to construct the discrete phase-type distribution. Please provide either n or the subintensity matrix.')
  }
  if(is.null(itons)){
    ph = t_total(n)

    T_table = ph$subint_mat
    alpha = ph$init_probs
    defect = 0

    P = solve(diag(nrow(T_table)) - 2/theta * T_table)
  }
  else if(itons <= 0 | itons > (n-1) | !is.numeric(itons)){
    stop('itons should be a number between 1 and n-1')

  }
  else if (is.matrix(subint_mat)) {
    if(!is.numeric(subint_mat) | nrow(subint_mat) != ncol(subint_mat)){
      stop('Subintensity matrix should be a square numerical matrix')
    }
    # rowsums in Subintensity matrix have to be non positive
    else if(sum(rowSums(subint_mat) > 0) != 0){
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
