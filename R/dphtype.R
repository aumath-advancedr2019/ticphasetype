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
#' So far only works for T_MRCA
#'
#' @usage rphtype(n, init_probs, subint_mat)
#'
#' @export

rphtype <- function(n_samples, n_OTU, type = "T_MRCA") {

  if (type == "T_MRCA") {
    # A: Calculate density function. Break when the number is very low
    init_probs = generate_init_row(n_OTU - 1)
    subint_mat = generate_subint_mat(n_OTU)

    x = 100000

    # I copied the dphtype function into here, because I need a contingent break in the loop.
    vec <- c()
    e <- matrix(rep(1, nrow(subint_mat)), nrow(subint_mat), 1)
    for (i in seq(2, x, 0.01)) {
      new_item =  -init_probs%*%expm(i*subint_mat)%*%subint_mat%*%e
      vec <- c(vec,new_item)
      if (i > 4 & new_item < 0.0000000001) {
        break # TODO: use the mean (calculate using PH) to know when to look for infinitesimal value. (instead of just `i>4`)
      }
    }
    print(length(vec))
    print(vec[1:100])


    # B: put the density function (stored in variable: `vec`) into a sampling function that accepts a weight parameter.
    sample(seq(2, x, 0.01)[1:length(vec)], n_samples, replace = T, prob = vec)
  }

}


#OTUs = 5

#present rph
#tibble(x = rphtype(100000, OTUs)) %>% ggplot(aes(x)) + geom_histogram(binwidth = 0.01) + labs(title = "10000 random samples from density with 5 OTUs")


#present dph
#tibble(x = seq(2, 10, 0.01),
#       p = dphtype(seq(2, 10, 0.01), generate_init_row(OTUs-1), generate_subint_mat(OTUs))) %>%
#  ggplot(aes(x, p)) + geom_line() + labs(title = "density for 5 OTUs")

#sample(seq(2, 100000, 0.01)) %>% hist





