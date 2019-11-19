
#' Simulate Matrices from Kingsman Coalescent
#'
#' Wrapper function that computes the matrices descrbining Kingsman coalescent ready for reward transformation
#'
#' @usage kingsman(n)
#'
#' @param n sample size (should be a positive integer)
#'
#' @return A `mult_phase_type` object containing the subintensity matrix, reward matrix and vector of initial probabilities
#'
#' @examples
#' kingsman(4)
#'
#' @export

kingsman <- function(n = NULL){

  if (is.null(n) | n <= 1 | !is.numeric(n) | n %% 1 != 0 | n %% 1 != 0) {
    stop('n should be a positive integer larger than 1')
  }

  else{
    matrixes = RateMAndStateSpace(n)
    # Rate Matrix
    if(n == 2){
      subint_mat = matrix(matrixes$RateM[1, 1])
    }
    else{
      subint_mat = matrixes$RateM[1:ncol(matrixes$RateM)-1, 1:ncol(matrixes$RateM)-1]
    }

    # Initial Distribution
    if(n == 2){
        init_probs = matrix(c(1))
      }
      else{
        init_probs = c(1, rep(0, nrow(subint_mat)-1))
      }
    }
  # Reward Matrix
  RewardM = matrixes$StSpM[1:(nrow(matrixes$StSpM)-1), 1:ncol(matrixes$StSpM)-1]

  value = list(subint_mat = subint_mat, RewardM = RewardM, init_probs = init_probs)
  attr(value, "class") <- c('mult_phase_type')
  return(value)
}

#' Discretize MPH object with a reward vector
#'
#' Reward transformation of continuous Phase Type distribution into Discrete Phase Type
#'
#' @usage RewTransform(mph_obj, rewards, theta)
#'
#' @param mph_obj Multivariate Phase Type object generated either from mult_phase_type or kingsman function
#' @param rewards vector of non negative numbers
#' @param theta mutation parameter (positive number)
#'
#' @return A `disc_phase_type` object containing subintensity matrix (P), vector of initial probabilities (alpha) and defect (probability of not entering any transient
#' state prior to absorption)
#'
#' @examples
#' RewTransform(kingsman(4), c(4, 2, 1, 0), 2)
#'
#' @export

RewTransform <- function(mph_obj, rewards = NULL, theta = NULL){
  if(is.null(rewards)){
    stop('rewards should be a 1D numerical vector')
  }
  else if(is.numeric(rewards) & sum(rewards < 0) != 0){
    stop('rewards has to be a vector of non negative numbers')
  }
  else{
    if(!is.numeric(rewards)) {
      if (class(rewards == 'matrix') & ncol != 1){
        stop('rewards should be a 1D numerical vector')
      }
      stop('rewards should be a 1D numerical vector')
    }
  }

  if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }
  ######### Computation of T*, alpha and defect ##########
  rew_transformed = rewardtransformparm(rewards, mph_obj$init_probs, mph_obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(subint_mat = P, init_probs = alpha, defect = defect)
  attr(value, 'class') <- c('disc_phase_type')
  value
}

#' Segregating Sites
#'
#' Generate subintensity matrix for the special case of segregating sites
#'
#' @usage segsites(n, theta)
#'
#' @param n sample size (positive integer)
#' @param theta mutation parameter (positive)
#'
#' @return A `disc_phase_type` object containing subintensity matrix (P), vector of initial probabilities (alpha) and defect (probability of not entering any transient
#' state prior to absorption)
#'
#' @examples
#' segsites(n = 4, theta = 2)
#'
#' @export

segsites <- function(n, theta=2){
  if (n<=1 | !is.numeric(n) | n %% 1 != 0 | n %% 1 != 0) {
    stop('n should be a positive integer larger than 1')
  }
  else if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }

  ph = t_total(n)

  T_table = ph$subint_mat
  alpha = ph$init_probs
  defect = 0

  P = solve(diag(nrow(T_table)) - 2/theta * T_table)

  value <- list(subint_mat = P, init_probs = alpha, defect = 0)
  attr(value, "class") <- c("disc_phase_type")
  return(value)
}
#' Itons (frequency counts)
#'
#' Generate subintensity matrix, vector of initial probabilites and defect for any frequency count
#'
#' @usage itons(mph_obj, itons, theta)
#'
#' @param mph_obj Multivariate Phase Type object generated either from mult_phase_type or kingsman function
#' @param itons Frequency count (positive integer)
#' @param theta mutation parameter (positive)
#'
#' @return A `disc_phase_type` object containing subintensity matrix (P), vector of initial probabilities (alpha) and defect (probability of not entering any transient
#' state prior to absorption)
#'
#' @examples
#' itons(kingsman(4), 2, theta = 2)
#'
#' @export

itons <- function(mph_obj, itons, theta) {
  if(!is.numeric(itons) | itons %% 1 != 0 | itons < 1){
    stop('itons should be a positive integer larger than 0')
  }
  else if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }
  else if(itons > ncol(mph_obj$RewardM)){
    stop('itons (frequency count) has to be lower than the sample size "n"')
  }
  if(length(mph_obj$subint_mat) == 1){
    reward = matrix(mph_obj$RewardM)
  }
  else{
    reward = mph_obj$RewardM[,itons]
  }
  ######### Computation of T*, alpha and defect ##########
  rew_transformed = rewardtransformparm(reward, mph_obj$init_probs, mph_obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(subint_mat = P, init_probs = alpha, defect = defect)
  attr(value, 'class') <- c('disc_phase_type')
  value
}
#' Tail Statistic
#'
#' Generate subintensity matrix, vector of initial probabilites and defect for combination of frequency counts higher than "k" (tail statistic)
#'
#' @usage tailstat(mph_obj, k, theta)
#'
#' @param mph_obj Multivariate Phase Type object generated either from mult_phase_type or kingsman function
#' @param k Minimum Frequency count (positive integer)
#' @param theta mutation parameter (positive)
#'
#' @return A `disc_phase_type` object containing subintensity matrix (P), vector of initial probabilities (alpha) and defect (probability of not entering any transient
#' state prior to absorption)
#'
#' @examples
#' tailstat(kingsman(15), k = 10, theta = 2)
#'
#' @export

tailstat <- function(mph_obj, k, theta) {
  if(!is.numeric(k) | k %% 1 != 0 | k < 1 | k >= ncol(mph_obj$RewardM)){
    stop('k should be a positive integer larger than 0 and smaller than the sample size')
  }
  else if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }

  if(nrow(mph_obj$subint_mat) == k){
    reward = matrix(mph_obj$RewardM[, ncol(mph_obj$RewardM)])
  }
  else{
    reward = rowSums(mph_obj$RewardM[, k:ncol(mph_obj$RewardM)])
  }

  ######### Computation of T*, alpha and defect ##########
  rew_transformed = rewardtransformparm(reward, mph_obj$init_probs, mph_obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(subint_mat = P, init_probs = alpha, defect = defect)
  attr(value, 'class') <- c('disc_phase_type')
  value
}

