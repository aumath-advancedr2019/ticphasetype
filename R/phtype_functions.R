#' Density, distribution function, quantile function and random generation for
#' phase-type distributions
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param obj an object of class \code{phase_type} or \code{disc_phase_type}.
#' @param n number of samples to draw.
#' @param granularity distance between numbers drawable
#'
#' @export


#' @import expm
#'


#' @describeIn dphtype
#'
#' Density function.
#'
#' @usage dphtype(x, obj)
#'
#' @export

dphtype <- function(x, obj){
  if (class(obj) == 'phase_type') {
    vec <- c()
    e <- matrix(rep(1,nrow(obj$subint_mat)), nrow(obj$subint_mat), 1)
    for (i in x) {
      vec <- c(vec, -obj$init_probs%*%expm(i*obj$subint_mat)%*%obj$subint_mat%*%e)
    }
    return(vec)
  } else if (class(obj) == 'disc_phase_type') {
    e = matrix(1, nrow = nrow(obj$subint_mat))
    t = e - obj$subint_mat %*% e
    dens_vec = c()
    for(i in x){
      dens_vec <- c(dens_vec, obj$init_probs %*% (obj$subint_mat %^% (i-1)) %*% t)
    }
    return(dens_vec)
  } else {
    stop("Please provide a 'phase_type' or a 'disc_phase_type' class.")
  }
}


#' @describeIn dphtype
#'
#' Quantile Function. Still not implemented for disc_phase_type.
#'
#' @usage qphtype(p, obj)
#'
#' @import stats
#'
#' @export

qphtype <- function(p, obj){
  if (class(obj) == 'phase_type') {
    vec <- c()
    inv <- function(y) uniroot(function(q) pphtype(q, obj)-y, c(0,20))$root[1]
    for (i in p) {
      vec <- c(vec, inv(i))
    }
    return(vec)
  } else if (class(obj) == 'disc_phase_type') {
    stop('Still not implemented for disc_phase_type.')
  } else {
    stop("Please provide a 'phase_type' or a 'disc_phase_type' class.")
  }
}


#' @describeIn dphtype
#'
#' Distribution function.
#'
#' @usage pphtype(q, obj)
#'
#' @export


pphtype <- function(q, obj){
  if (class(obj) == 'phase_type') {
    vec <- c()
    e <- matrix(rep(1,nrow(obj$subint_mat)), nrow(obj$subint_mat), 1)
    for (i in q) {
      vec <- c(vec, 1-obj$init_probs%*%expm(i*obj$subint_mat)%*%e)
    }
    return(vec)
  } else if (class(obj) == 'disc_phase_type') {
    e = matrix(1, nrow = nrow(obj$subint_mat))
    prob_vec = c()
    for(i in q){
      prob_vec <- c(prob_vec, 1 - obj$init_probs %*% (obj$subint_mat %^% i) %*% e)
    }
    return(prob_vec)
  } else {
    stop("Please provide a 'phase_type' or a 'disc_phase_type' class.")
  }
}


#' @describeIn dphtype
#'
#' Random number generator. Still not implemented for disc_phase_type.
#'
#' @usage rphtype(n, obj, granularity = 0.01)
#'
#' @export


rphtype <- function(n, obj, granularity = 0.01){
  if (class(obj) == 'phase_type') {
    # A: Calculate density function. Break when the number is very low


    x = 100000

    # I copied the dphtype function into here, because I need a contingent break in the loop.
    vec <- c()
    e <- matrix(rep(1, nrow(obj$subint_mat)), nrow(obj$subint_mat), 1)
    for (i in seq(0, x, granularity)) {
      new_item =  -obj$init_probs%*%expm(i*obj$subint_mat)%*%obj$subint_mat%*%e
      vec <- c(vec, new_item)
      if (i > 4 & new_item < 0.0000000001) { # Only if you sample more than a billion, will you see a bias induced by the `break``
        break # TODO: use the mean (calculate using PH) to know when to look for infinitesimal value. (instead of just `i>4`)
      }
    }
    #print(length(vec))
    #print(vec[1:100])


    # B: put the density function (stored in variable: `vec`) into a sampling function that accepts a weight parameter.
    return(sample(seq(0, x, 0.01)[1:length(vec)], n, replace = T, prob = vec))
  } else if (class(obj) == 'disc_phase_type') {
    stop('Still not implemented for disc_phase_type.')
  } else {
    stop("Please provide a 'phase_type' or a 'disc_phase_type' class.")
  }
}












