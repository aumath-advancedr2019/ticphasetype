#' Time until the most recent common ancestor
#'
#' The function creates a continuous phase-type representation of the time until the most recent
#' common ancestor for a specified number of sequences and following Kingsman's coalescent model.
#'
#'
#' @param n positive integer larger than 1.
#'
#' @usage t_mrca(n)
#'
#' @examples
#' ph_mrca <- t_mrca(4)
#' summary(ph_mrca)
#'
#' @export


t_mrca <- function(n) {
  if (n<=1 | !is.numeric(n)) {
    stop('n must be an integer larger than 1.')
  }
  subint_mat = matrix(c(0), nrow = n-1, ncol = n-1)
  for (i in 1:n-1) {
    subint_mat[i, i] = - (n - i+1)*(n - i)/2
    if (i < n-1) {
      subint_mat[i, i+1] = -subint_mat[i, i]
    }
  }
  init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
  return(cont_phase_type(subint_mat = subint_mat, init_probs = init_probs))
}
