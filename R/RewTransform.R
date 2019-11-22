#' Discretize MPH object with a reward vector
#'
#' Reward transformation of continuous phase-type distribution into a discrete phase-type distribution.
#'
#' @usage RewTransform(mph_obj, rewards, theta)
#'
#' @param mph_obj multivariate phase-type object generated either from mult_phase_type or kingsman function
#' @param rewards vector of non negative numbers
#' @param theta mutation parameter (positive number)
#'
#' @return A `disc_phase_type` object containing subintensity matrix (P), vector of initial probabilities (alpha) and defect (probability of not entering any transient
#' state prior to absorption)


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
  rew_transformed = rewardtransformation(rewards, mph_obj$init_probs, mph_obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(subint_mat = P, init_probs = alpha, defect = defect)
  attr(value, 'class') <- c('disc_phase_type')
  value
}
