#' \code{disc_phase_type} class
#'
#' Description of the class \code{disc_phase_type}, which represents discrete phase-type
#' distributions.
#'
#' \code{disc_phase_type} is the generator function for the discrete phase-type distribution class
#' of the same name, which inherits from \code{list}. This function can easily produce phase-type
#' representations of common statistics in genomics, such as the total number of segregating sites,
#' the number of i-tons (singletons, doubletons, etc.) and the tail statistic.
#'
#' The discrete phase-type distribution can be generated in two different ways:
#' \itemize{
#'   \item By specifying whether the user wants a  discrete phase-type representation
#'   of the total number of segregating sites (default), the number of i-tons, or the tail
#'   statistic for a certain number of sequences \code{n}.
#'   \item By supplying a user-defined sub-intensity matrix, with optional initial
#'   probabilities.
#' }
#'
#' See examples for further explanation on the usage.
#'
#' @param object an object of class \code{disc_phase_type}.
#' @param n integer larger than 1, or \code{NULL} (default).
#' @param itons integer between 1 and n-1, or \code{NULL} (default).
#' @param init_probs vector, a one-row matrix or \code{NULL} (default).
#' @param subint_mat matrix or \code{NULL} (default).
#' @param theta numeric.
#' @param tail_stat logical.
#'
#' @usage disc_phase_type(n = NULL, itons = NULL,
#'             subint_mat = NULL, init_probs = NULL,
#'             theta = 2, tail_stat = F)
#'
#' @examples
#' # Total number of segregating sites
#' disc_phase_type(4)
#' disc_phase_type(5, theta=0.5)
#'
#' # Below examples are not implemented:
#' # Number of singletons
#' disc_ph_example <- disc_phase_type(4, itons=1)
#' mean(disc_ph_example)
#' var(disc_ph_example)
#' summary(disc_ph_example)
#'
#' # Tail statistic
#' disc_ph_tail <- disc_phase_type(4, itons=2, tail_stat=TRUE)
#' summary(disc_ph_tail)
#'
#' @export

disc_phase_type = function(n = NULL, itons = NULL, theta = 2, tail_stat = F, subint_mat = NULL, init_probs = NULL){

  if (is.null(n)) {
    if (is.null(subint_mat)) {
      stop('Unable to construct the discrete phase-type distribution. Please provide either n or the subintensity matrix.')
    } else if (is.matrix(subint_mat)) {
      if (is.null(init_probs)) {
        init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
        warning('The initial probability vector is automatically generated.')
      } else if ((is.vector(init_probs) & is.atomic(init_probs)) | is.matrix(init_probs)) {
        if (nrow(subint_mat) == length(init_probs)) {
          init_probs <- matrix(init_probs, nrow = 1)
        } else {
          stop('The length of the initial probabilities does not match the size of the subintensity matrix.')
        }
      } else {
        stop('The initial probabilities must be a a matrix with one row or a vector.')
      }
    } else {
      stop('The subintensity matrix must be a matrix.')
    }
    value = list(subint_mat = subint_mat, init_probs = init_probs, defect = 1-sum(init_probs))
    attr(value, "class") <- "phase_type"
    return(value)
  }


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
  attr(value, 'class') <- 'disc_phase_type'
  value
}


#' @describeIn disc_phase_type
#'
#' mean of the discrete phase-type distribution.
#'
#' @usage ## S3 method for class 'disc_phase_type'
#' mean(object)
#'
#' @export

mean.disc_phase_type <- function(obj) {
  mean <- sum(obj$init_probs%*%solve(diag(nrow = nrow(obj$subint_mat))-obj$subint_mat))
  as.numeric(mean+obj$defect)
}

#' @export

var <- function(x, ...) {
  UseMethod('var', x)
}

#' @describeIn disc_phase_type
#'
#' variance of the discrete phase-type distribution.
#'
#' @usage ## S3 method for class 'disc_phase_type'
#' var(object)
#'
#' @export

var.disc_phase_type <- function(obj) {
  variance <- sum(2*obj$init_probs%*%obj$subint_mat%*%solve((diag(nrow = nrow(obj$subint_mat))-obj$subint_mat)%^%2)) +
    mean(obj) -
    mean(obj)^2
  as.numeric(variance)
}

#' @describeIn disc_phase_type
#'
#' summary of the discrete phase-type distribution.
#'
#' @usage ## S3 method for class 'disc_phase_type'
#' summary(object)
#'
#' @export

summary.disc_phase_type <- function(obj) {
  cat('\nSubintensity matrix:\n')
  print(obj$subint_mat)
  cat('\nInitial probabilities:\n')
  print(obj$init_probs)
  cat('\nDefect:\n')
  print(obj$defect)
  cat('\nMean: ', mean(obj), '\n', sep = '')
  cat('\nVariance: ', var(obj), '\n\n', sep = '')
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

