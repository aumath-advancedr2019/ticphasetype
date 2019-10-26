#' Density, distribution function, quantile function and random generation for
#' the continuous phase-type distribution
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param init_probs vector of initial probabilities.
#' @param subint_mat subintensity matrix.
#'
#' @export


#' @describeIn dphtype
#'
#' Density function.
#'
#' @usage dphtype(x, init_probs, subint_mat)
#'
#' @export

dphtype <- function(x, init_probs, subint_mat) {
  vec <- c()
  e <- matrix(rep(1,nrow(subint_mat)), nrow(subint_mat), 1)
  for (i in x) {
    vec <- c(vec, -init_probs%*%expm(i*subint_mat)%*%subint_mat%*%e)
  }
  vec
}

#' @describeIn dphtype
#'
#' Quantile Function.
#'
#' @usage qphtype(q, init_probs, subint_mat)
#'
#' @export

qphtype <- function(p, init_probs, subint_mat) {
  vec <- c()
  inv <- function(y) uniroot(function(q) pphtype(q, init_probs, subint_mat)-y, c(0,20))$root[1]
  for (i in p) {
    vec <- c(vec, inv(i))
  }
  vec
}

#' @describeIn dphtype
#'
#' Distribution function.
#'
#' @usage pphtype(p, init_probs, subint_mat)
#'
#' @export

pphtype <- function(q, init_probs, subint_mat) {
  vec <- c()
  e <- matrix(rep(1,nrow(subint_mat)), nrow(subint_mat), 1)
  for (i in q) {
    vec <- c(vec, 1-init_probs%*%expm(i*subint_mat)%*%e)
  }
  vec
}

#' @describeIn dphtype
#'
#' Random number generator.
#'
#' @usage rphtype(n, init_probs, subint_mat)
#'
#' @export

rphtype <- function(n, init_probs, subint_mat) {

}





