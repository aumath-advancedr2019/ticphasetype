#' Site Frequency Spectrum
#'
#' @description Returns a list of expected counts of itons up to n-1 or a plot if specified. Works with different values of `theta` (population mutation parameter)
#' as well as with custom vector of initial probabilities.
#'
#' @usage sfs(n, theta = 2, init_probs = NULL, plot = FALSE)
#'
#' @param n sample size
#' @param theta population mutation parameter
#' @param init_probs vector of initial probabilities. There is no need to supply entire vector if most of its elements are zero. The function computes the rest automatically
#' @param plot logical
#'
#' @examples sfs(n = 5)
#' sfs(5, theta = 3)
#' sfs(10, 2, init_probs = c(0.8, 0.2))
#'
#' @import graphics
#'
#' @export
sfs <- function(n, theta = 2, init_probs = NULL, plot = FALSE){
  E_ksi = rep(NA, n-1)
  Var_ksi = rep(NA, n-1)

  states = nrow(RateMAndStateSpace(n)$StSpM)-1
  if(!is.null(init_probs)){
    if(length(init_probs) != states){
      message('Zeros added to init_probs vector to make it the right size')
      init_probs = c(init_probs, rep(0, states - length(init_probs)))
    }
  }
  # calculate E[ksi_i] and Var[ksi_i] for every i < n
  for (iton in 1:(n-1)){
    dph = disc_phase_type(n, itons = iton, theta = theta, init_probs = init_probs)
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
