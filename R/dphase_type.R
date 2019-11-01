
#' \code{dphase_type} class
#'
#' Description of the class \code{phase_type}, which represents discrete phase-type
#' distributions.
#'
#' @param object an object of class \code{dphase_type}.
#' @param type \code{'T_MRCA'}, \code{'T_Total'} or \code{NULL} (default).
#' @param n integer larger than 1.
#' @param itons integer between 1 and n-1, or \code{NULL} (default).
#' @param init_probs vector, a one-row matrix or \code{NULL} (default).
#' @param theta numeric.
#'
#' @usage dphase_type(n, init_probs = NULL, itons = NULL, theta = 2, moment = F, tail_stat = F)
#'
#'
#' @export

dphase_type = function(n, init_probs = NULL, itons = NULL, theta = 2, moment = F, tail_stat = F){

  # All Segregating Sites - This computes the P matrix faster than setting rewards to sum(ri)
  if (n<=1 | !is.numeric(n)) {
    stop('n should be a positive integer larger than 1')
  }

  if(is.null(itons)){
    ph = phase_type('T_Total', n = n)

    T_table = ph$subint_mat
    alpha = ph$init_probs
    defect = 0

    P = solve(diag(nrow(T_table)) - 2/theta * T_table)
  }
  else if(itons <= 0 | itons > (n-1) | !is.numeric(itons)){
    stop('itons should be a number between 1 and n-1')
  }
  # Special Case (eg. singletons or tail statistic)
  else{

    ######## Step1: Preparation of Rate matrix (T), initial Distribution (pi) and reward vector ########
    matrixes = RateMAndStateSpace(n)
    # Rate Matrix
    if(n == 2){
      T_table = matrix(matrixes$RateM[1, 1])
    }
    else{
      T_table = matrixes$RateM[1:ncol(matrixes$RateM)-1, 1:ncol(matrixes$RateM)-1]
    }

    # Initial Distribution
    if(is.null(init_probs)){
      if(n == 2){
        init_probs = matrix(c(1))
      }
      else{
        init_probs = c(1, rep(0, nrow(T_table)-1))
      }
    }

    ####### REWARDS #######
    # Tail Statistic
    if(tail_stat){
      if(n == (itons + 1)){
        reward = matrix(matrixes$StSpM[1:(nrow(matrixes$StSpM)-1), ncol(matrixes$StSpM)-1])
      }
      else{
        reward = apply(matrixes$StSpM[1:(nrow(matrixes$StSpM)-1), itons:(ncol(matrixes$StSpM)-1)], 1, sum)
      }
    }

    # I-tons
    else{
      if(n == 2){
        reward = matrix(matrixes$StSpM[1,1])
      }
      else{
        reward = matrixes$StSpM[1:nrow(matrixes$StSpM)-1,][,itons]
      }
    }

    ######### Step2: Computation of T*, alpha and defect ##########
    rew_transformed = rewardtransformparm(reward, init_probs, T_table)
    alpha = rew_transformed$init_probs
    T_star = rew_transformed$subint_mat
    defect = rew_transformed$defect

    ########## Step3: Computation of P and p (transformation to DPH) ##########
    P = solve(diag(nrow(T_star)) - 2/theta * T_star)
  }

  if(is.null(itons)){
    value = list(subint_mat = P, init_probs = alpha, defect = defect)
  }
  else{
    value = list(subint_mat = P, init_probs = alpha, defect = defect)
  }
  attr(value, 'class') <- 'dphase_type'
  value
}


#' @describeIn dphase_type
#'
#' mean of the discrete phase-type distribution.
#'
#' @usage ## S3 method for class 'dphase_type'
#' mean(object)
#'
#' @export

mean.dphase_type <- function(obj) {
  mean <- sum(obj$init_probs%*%solve(diag(nrow = nrow(obj$subint_mat))-obj$subint_mat))
  as.numeric(mean+obj$defect)
}

#' @export

var <- function(x, ...) {
  UseMethod('var', x)
}

#' @describeIn dphase_type
#'
#' variance of the discrete phase-type distribution.
#'
#' @usage ## S3 method for class 'dphase_type'
#' var(object)
#'
#' @export

var.dphase_type <- function(obj) {
  variance <- sum(2*obj$init_probs%*%obj$subint_mat%*%solve((diag(nrow = nrow(obj$subint_mat))-obj$subint_mat)%^%2)) +
    mean(obj) -
    mean(obj)^2
  as.numeric(variance)
}

#' @describeIn dphase_type
#'
#' summary of the discrete phase-type distribution.
#'
#' @usage ## S3 method for class 'dphase_type'
#' summary(object)
#'
#' @export

summary.dphase_type <- function(obj){
  cat('Mean:\n',mean(obj))
  cat('\nVariance:\n',var(obj))
}

sfs <- function(n_vec, theta = 2){
  sfs_df <- tibble()
  # iterate through sample sizes given in n_vec
  for (i in 1:length(n_vec)){

    # calculate E[ksi_i] and Var[ksi_i] for every i < n
    for (iton in 1:(n_vec[i]-1)){
      dph = discph(n_vec[i], itons = iton, moment = T, theta = theta)
      varksi = var(dph)
      ksi = dmoments(dph, 1) - 1 + dph$defect
      sfs_df <- bind_rows(sfs_df, tibble(n = n_vec[i], itons = iton, E_ksi = ksi, Var_ksi = varksi))
    }
  }
  sfs_df
}

#----------------------------------------------------------------------------------------------
#' dsegsites
#'
#' Distribution function for segregating sites
#'
#' @usage pdphtype(k, alpha, T_star, theta)
#'
#' @export

dsegsites <- function(x, obj){
  p = matrix(1, nrow = nrow(obj$subint_mat)) - P %*% matrix(1, nrow = nrow(obj$subint_mat))
  dens_vec = c()
  for (i in x){
    dens_vec <- c(dens_vec, obj$init_probs %*% (obj$subint_mat%^%(i)) %*% p)
  }
  dens_vec
}

