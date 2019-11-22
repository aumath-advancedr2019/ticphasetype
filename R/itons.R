#' I-tons (frequency counts)
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
  rew_transformed = rewardtransformation(reward, mph_obj$init_probs, mph_obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  disc_phase_type(P, alpha)
}
