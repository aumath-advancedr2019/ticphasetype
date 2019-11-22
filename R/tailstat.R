



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
  rew_transformed = rewardtransformation(reward, mph_obj$init_probs, mph_obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  disc_phase_type(P, alpha)
}

