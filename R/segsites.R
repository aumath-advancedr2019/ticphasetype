#' Segregating sites
#'
#' Generate subintensity matrix for the special case of segregating sites
#'
#' @usage segsites(n, theta)
#'
#' @param n sample size (positive integer)
#' @param theta mutation parameter (positive)
#'
#' @return A `disc_phase_type` object containing subintensity matrix (P), vector of initial probabilities (alpha) and defect (probability of not entering any transient
#' state prior to absorption)
#'
#' @examples
#' segsites(n = 4, theta = 2)
#'
#' @export

segsites <- function(n, theta=2){
  if (n<=1 | !is.numeric(n) | n %% 1 != 0 | n %% 1 != 0) {
    stop('n should be a positive integer larger than 1')
  }
  else if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }

  ph = t_total(n)

  T_table = ph$subint_mat
  alpha = ph$init_probs
  defect = 0

  P = solve(diag(nrow(T_table)) - 2/theta * T_table)

  value <- list(subint_mat = P, init_probs = alpha, defect = 0)
  attr(value, "class") <- c("disc_phase_type")
  return(value)
}
