#' Site Frequency Spectrum
#'
#' This function returns a list of expected counts of itons up to n-1, or a plot if specified.
#' It builds the site frequency spectrum using Kingman's coalescent. It works for different
#' values of `theta` (population mutation parameter).
#'
#' @usage sfs(n, theta = 2, plot = FALSE)
#'
#' @param n sample size
#' @param theta population mutation parameter
#' @param plot logical
#'
#' @details
#'
#' By default, a list of expected values for each i-ton and its variance is returned.
#' If \code{plot=TRUE}, then a plot is returned.
#'
#' @examples
#' sfs(n = 5, theta = 3)
#'
#' @import graphics
#'
#' @export


sfs <- function(n, theta, plot = FALSE){
  E_ksi = rep(NA, n-1)
  Var_ksi = rep(NA, n-1)

  coalescent = kingsman(n)
  # calculate E[ksi_i] and Var[ksi_i] for every i < n
  for (iton in 1:(n-1)){
    dph = itons(coalescent, iton, theta)
    E_ksi[iton] = mean(dph) - 1
    Var_ksi[iton] = var(dph)
  }

  if (plot == TRUE) {
    data = list(E_ksi = E_ksi, Var_ksi = Var_ksi)
    barplot(data$E_ksi, names.arg = 1:(n-1),
            main = paste0('SFS for n=',n,' and theta=',theta),
            xlab = 'i-tons', ylab = 'Counts')
  } else {
    return(data = list(E_ksi = E_ksi, Var_ksi = Var_ksi))
  }


}
