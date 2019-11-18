#' Total tree length
#'
#' This function creates a continuous phase-type representation of the total tree length
#' for a specified number of sequences and following Kingsman's coalescent model.
#'
#' @param n positive integer larger than 1.
#'
#' @usage t_total(n)
#'
#' @examples
#' ph_total <- t_total(4)
#' summary(ph_total)
#'
#' @export

t_total <- function(n) {
  if (n<=1 | !is.numeric(n)) {
    stop('n must be an integer larger than 1.')
  }
  subint_mat = matrix(c(0), nrow = n-1, ncol = n-1)
  for (i in 1:n-1) {
    subint_mat[i, i] = - 0.5 * (n - i)
    if (i < n-1) {
      subint_mat[i, i+1] = -subint_mat[i, i]
    }
  }
  init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
  return(cont_phase_type(subint_mat = subint_mat, init_probs = init_probs))
}
