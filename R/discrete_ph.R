dphase_type <- function(n, subrate_mat = NULL, init_probs = NULL, rewards = NULL, itons = NULL, theta = 2){

  ############## Step1: Preparation of Rate matrix (T), initial probabilities (pi) and reward vector ##############

  #### Rate Matrix ###
  if(is.null(subrate_mat)){
    matrixes = RateMAndStateSpace(n)

    if(n == 2) {subrate_mat = matrix(matrixes$RateM[1, 1])}
    else{subrate_mat = matrixes$RateM[1:ncol(matrixes$RateM)-1, 1:ncol(matrixes$RateM)-1]}
  }

  #### Initial Distribution (init_probs) ####
  # if nothing is supplied in the function argument, it creates a vector with first entry
  # being 1 and all the others 0
  if(is.null(init_probs)){
    if(n == 2){init_probs = matrix(c(1))}
    else{init_probs = c(1, rep(0, nrow(subrate_mat)-1))}
  }

  ### REWARDS ###
  if(is.null(rewards)){
    # Specifying if we are considering all segregating sites or something more specific
    if(is.null(itons)){ # means all (singletons + doubletons + ...)
      if(n == 2){reward = matrix(sum(matrixes$StSpM[1,]))}
      else{reward = apply(matrixes$StSpM[1:nrow(matrixes$StSpM)-1, ], 1, sum)}
    }
    else{
      if(itons < n){
        if(n == 2){reward = matrix(matrixes$StSpM[1,1])}
        else{reward = matrixes$StSpM[1:nrow(matrixes$StSpM)-1,][,itons]}
      }
      else{return(0)}
    }
  }

  ############## Step2: Computation of T*, alpha and defect ##############
  rew_transformed = rewardtransformparm(reward, init_probs, subrate_mat)
  alpha = rew_transformed$newinitprob
  T_star = rew_transformed$newsubintensitymatrix

  ############## Step3: Computation of subtransition matrix P ##############
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(P = P, T_star = T_star, alpha = alpha, defect = rew_transformed$defect)
  attr(value, "class") <- "dphase_type"
  value
}


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

#' @import matrixcalc

dmoments <- function(stm, ip, m){
  e = matrix(1, nrow = nrow(stm))
  moment = factorial(m) * ip %*% (stm %^% (m-1)) %*% matrixcalc::matrix.power(diag(nrow = nrow(stm)) - stm, -m) %*% e
  #e = matrix(1, nrow = nrow(phase_type$subint_mat))
  #moment = factorial(m) * phase_type$init_probs %*% (phase_type$subint_mat %^% (m-1)) %*% matrixcalc::matrix.power(diag(nrow = nrow(phase_type$subint_mat)) - phase_type$subint_mat, -m) %*% e
  moment[1]
}

dvariance <- function(phase_type){
  m1 = dmoments(phase_type, 1)
  m2 = dmoments(phase_type, 2)
  m2 - m1^2 + m1
}
