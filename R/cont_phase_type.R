#' \code{cont_phase_type} class
#'
#' Description of the class \code{cont_phase_type}, which represents continuous phase-type
#' distributions.
#'
#' \code{cont_phase_type} is the generator function for the continuous phase-type distribution class
#' of the same name, which inherits from \code{list}. It can be generated in two different ways:
#' \itemize{
#'   \item By specifying whether the user wants a phase-type representation of the
#'   total time until the most recent common ancestor (MRCA) or the total tree
#'   length, for a certain number of sequences \code{n}.
#'   \item By supplying a user-defined sub-intensity matrix, with optional initial
#'   probabilities.
#' }
#'
#' @param type \code{'T_MRCA'}, \code{'T_Total'} or \code{NULL} (default).
#' @param n integer larger than 1 or \code{NULL} (default).
#' @param subint_mat matrix or \code{NULL} (default).
#' @param init_probs vector, a one-row matrix or \code{NULL} (default).
#'
#' @usage cont_phase_type(type = NULL, n = NULL,
#'            subint_mat = NULL, init_probs = NULL)
#'
#' @examples
#' #library(qpdf) # I added this for devtools::check() to pass.
#' # Time until the MRCA
#' cont_phase_type('T_MRCA', 4)
#'
#' # Total tree length
#' ph_example <- cont_phase_type('T_Total', 6)
#' mean(ph_example)
#' var(ph_example)
#' summary(ph_example)
#'
#' # User-specified phase-type distribution
#' subint_example <- matrix(runif(16), ncol=4)
#' ph_user <- cont_phase_type(subint_mat = subint_example)
#' summary(ph_user)
#'
#' @export



cont_phase_type <- function(type = NULL, n = NULL, subint_mat = NULL, init_probs = NULL) {
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
  value <- list(subint_mat = subint_mat,
                init_probs = init_probs,
                defect = 1-sum(init_probs))
  attr(value, "class") <- "cont_phase_type"
  value
}


moment_ph <- function(obj, m) {
  e <- matrix(rep(1,nrow(obj$subint_mat)), nrow(obj$subint_mat), 1)
  inv <- solve(obj$subint_mat%^%m)
  as.numeric((-1)**m*factorial(m)*obj$init_probs%*%inv%*%e)
}



#' @export

mean.cont_phase_type <- function(x, ...) {
  moment_ph(x, 1)
}

#' Variance of cont_phase_type distributions
#'
#' It calculates the variance of continuous and discrete phase-type distributions,
#' represented by the \code{cont_phase_type} and the \code{disc_phase_type} classes
#' respectively
#'
#' @param obj a cont_phase_type or disc_phase_type object.
#'
#' @export

var <- function(obj) {
  UseMethod('var', obj)
}

#' @export

var.cont_phase_type <- function(obj) {
  moment_ph(obj, 2)-moment_ph(obj, 1)**2
}


#' @export

summary.cont_phase_type <- function(object, ...) {
  cat('\nSubintensity matrix:\n')
  print(object$subint_mat)
  cat('\nInitial probabilities:\n')
  print(object$init_probs)
  cat('\nDefect:\n')
  print(object$defect)
  cat('\nMean: ', mean(object), '\n', sep = '')
  cat('\nVariance: ', var(object), '\n\n', sep = '')
}







