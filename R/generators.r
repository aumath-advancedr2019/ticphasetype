#' Helper functions to generate inputs for the phasetype-related functions.
#'
#' @param n number of samples.
#' @param type string denoting the type of PH-related subintensity matrix wanted: T_MRCA, S_total
#'
#' @export






#' @describeIn generate_subint_mat
#'
#' Generates a subintensity matrix from a specified number of samples.
#' Currently only supports generating the subintensity matrix used in calculation
#' of the time to the most recente common ancestor (MRCA)
#'
#'
#' @usage rphtype(n, init_probs, subint_mat)
#'
#' @export

generate_subint_mat <- function(n, type = "T_MRCA") {
  if (type == "T_MRCA") {
    T = matrix(c(0), nrow = n-1, ncol = n-1)
    for (i in 1:n-1) {
      #T[i, i] = - choose(n-i+1, 2)
      T[i, i] = - (n - i+1)*(n - i)/2 # equivalent to the above

      if (i < n-1) {
        T[i, i+1] = -T[i, i]
      }
    }
  }

  T
}

generate_init_row <- function(n) {
  c(1, rep(0, n-1))
}
