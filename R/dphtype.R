#' Density, distribution function, quantile function and random generation for
#' the continuous phase-type distribution
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param phase_type an object of class \code{phase_type}.
#' @param n_samples number of samples to draw.
#' @param n_OTU number of OTUs to consider
#' @param granularity distance between numbers drawable
#' @param type type of calculation ("T_MRCA" and later S_total)
#'
#' @export


#' @import expm
#'



#' @describeIn dphtype
#'
#' Density function.
#'
#' @usage dphtype(x, phase_type)
#'
#' @export

dphtype <- function(x, phase_type) {
  vec <- c()
  e <- Matrix(rep(1,nrow(phase_type$subint_mat)), nrow(phase_type$subint_mat), 1, sparse = T)
  for (i in x) {
    vec <- c(vec, (-phase_type$init_probs%*%expm.Higham08(i*phase_type$subint_mat, balancing = F)%*%phase_type$subint_mat%*%e)[1,1])
  }
  vec
}

#' @describeIn dphtype
#'
#' Quantile Function.
#'
#' @usage qphtype(q, phase_type)
#'
#' @export

qphtype <- function(p, phase_type) {
  vec <- c()
  inv <- function(y) uniroot(function(q) pphtype(q, phase_type)-y, c(0,20))$root[1]
  for (i in p) {
    vec <- c(vec, inv(i))
  }
  vec
}

#' @describeIn dphtype
#'
#' Distribution function.
#'
#' @usage pphtype(p, phase_type)
#'
#' @export

pphtype <- function(q, phase_type) {
  vec <- c()
  e <- Matrix(rep(1,nrow(phase_type$subint_mat)), nrow(phase_type$subint_mat), 1, sparse = T)
  for (i in q) {
    vec <- c(vec, (1-phase_type$init_probs%*%expm.Higham08(i*phase_type$subint_mat, balancing = F)%*%e)[1,1])
  }
  vec
}

#' @describeIn dphtype
#'
#' Random number generator.
#' So far only works for T_MRCA
#'
#' @usage rphtype(n_samples, phase_type, granularity = 0.01)
#'
#' @export

rphtype <- function(n_samples, phase_type, granularity = 0.01) {


  # A: Calculate density function. Break when the number is very low


  x = 100000

  # I copied the dphtype function into here, because I need a contingent break in the loop.
  vec <- c()
  e <- Matrix(rep(1, nrow(phase_type$subint_mat)), nrow(phase_type$subint_mat), 1, sparse = T)
  for (i in seq(0, x, granularity)) {
    new_item =  (-phase_type$init_probs%*%expm.Higham08(i*phase_type$subint_mat, balancing = F)%*%phase_type$subint_mat%*%e)[1,1]
    vec <- c(vec, new_item)
    if (i > 4 & new_item < 0.0000000001) { # Only if you sample more than a billion, will you see a bias induced by the `break``
      break # TODO: use the mean (calculate using PH) to know when to look for infinitesimal value. (instead of just `i>4`)
    }
  }
  #print(length(vec))
  #print(vec[1:100])


  # B: put the density function (stored in variable: `vec`) into a sampling function that accepts a weight parameter.
  sample(seq(0, x, 0.01)[1:length(vec)], n_samples, replace = T, prob = vec)


}








