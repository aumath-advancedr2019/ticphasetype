#' Density, distribution function, quantile function and random generation for
#' the discrete phase-type distribution
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
# #' @param phase_type an object of class \code{phase_type}.
#' @param n_samples number of samples to draw.
#' @param n_OTU number of OTUs to consider
#' @param granularity distance between numbers drawable
#' @param type type of calculation ("T_MRCA" and later S_total)
#'
#' @export

#' @import expm
#'

#' @describeIn discrete_ph
#'
#' Density function.
#'
#' @usage ddphtype(x, subint_mat, init_probs)
#'
#' @export

ddphtype <- function(x, dphase_type){
  e = matrix(1, nrow = nrow(phase_type$subint_mat))
  t = e - phase_type$subint_mat %*% e
  dens_vec = c()
  for(i in x){
    dens_vec <- c(dens_vec, phase_type$init_probs %*% (phase_type$subint_mat %^% (i-1)) %*% t)
  }
  dens_vec
}

#' @describeIn discrete_ph
#'
#' Distribution function.
#' @usage pdphtype(x, subint_mat, init_probs)
#'
#' @export

pdphtype <- function(q, phase_type){
  e = matrix(1, nrow = nrow(phase_type$subint_mat))
  prob_vec = c()
  for(i in q){
    prob_vec <- c(prob_vec, 1 - phase_type$init_probs %*% (phase_type$subint_mat %^% i) %*% e)
  }
  prob_vec
}

