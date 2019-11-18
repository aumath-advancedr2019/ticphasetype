cont_cont_phase_type(subint_mat, init_probs)
disc_phase_type(subint_mat, init_probs)
mult_phase_type(subint_mat, init_probs, rewards)

kingsman(n)
itons(mult_cont_phase_type, itons, theta)
tailstat(mult_cont_phase_type, tail, theta)
segsites(n, theta)
t_mrca(n)
t_total(n)

# obj: cont_cont_phase_type --> reward-transformed disc_phase_type
RewTransform(obj, reward)
RateMAndStateSpace(n)




#' @export

# This function creates a discrete phase-type distribution of Kingsman's coalescent
# together with the reward matrix for calculating the site-frequency spectrum and
# related statistics. The user can also specify a different coalescent model by
# supplying a sub-intensity matrix and an optional initial probabilities vector.
discrete_ph <- function(n = NULL, subint_mat = NULL, init_probs = NULL){

  if (is.null(n)) {
    if (is.null(subint_mat)) {
      stop('Unable to construct the discrete phase-type distribution. Please provide either n or the subintensity matrix.')
    }
    else if (is.matrix(subint_mat)) {
      if(!is.numeric(subint_mat) | nrow(subint_mat) != ncol(subint_mat)){
        stop('Subintensity matrix should be a square numerical matrix')
      }
      # i am not sure if the condition that values off diagonal have to be non negative and the ones on diagonal have to be
      # negative has to be neccessarily true.

      # rowsums in Subintensity matrix have to be non positive
      else if(sum(rowSums(subint_mat) > 0) != 0){
        stop('The rowsums in subintensity matrix have to be non-positive')
      }

      if (is.null(init_probs)) {
        init_probs <- matrix(c(1, rep(0, nrow(subint_mat) - 1)), 1, nrow(subint_mat))
        warning('The initial probability vector is automatically generated.')
      }
      else if ((is.vector(init_probs) & is.atomic(init_probs)) | is.matrix(init_probs)) {
        if (nrow(subint_mat) == length(init_probs)) {
          init_probs <- matrix(init_probs, nrow = 1)
        }
        else {
          stop('The length of the initial probabilities does not match the size of the subintensity matrix.')
        }
      }
      else {
        stop('The initial probabilities must be a a matrix with one row or a vector.')
      }
    }
    else {
      stop('The subintensity matrix must be a matrix.')
    }
    value = list(subint_mat = subint_mat, init_probs = init_probs, defect = 1-sum(init_probs))
    attr(value, "class") <- "disc_phase_type"
    return(value)
  }

  else if (n<=1 | !is.numeric(n) | n %% 1 != 0 | n %% 1 != 0) {
    stop('n should be a positive integer larger than 1')
  }
  else if(!is.null(init_probs)){
    if(!is.numeric(init_probs)) {
      if (class(init_probs == 'matrix') & ncol != 1){
        stop('init_probs should be a 1D numerical vector')
      }
      stop('init_probs should be a 1D numerical vector')
    }
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
    if(is.null(init_probs)){
      if(n == 2){
        init_probs = matrix(c(1))
      }
      else{
        init_probs = c(1, rep(0, nrow(subint_mat)-1))
      }
    }
  }
  # Reward Matrix
  RewardM = matrixes$StSpM[1:(nrow(matrixes$StSpM)-1), 1:ncol(matrixes$StSpM)-1]

  value = list(subint_mat = subint_mat, RewardM = RewardM, init_probs = init_probs, defect = 1-sum(init_probs))
  attr(value, "class") <- c('cont_phase_type', 'rewards')
  return(value)
}

# We need to have this function, otherwise there is no reason for having an `subint_mat` argument in the disc-ph function above
RewTransform <- function(obj, rewards = NULL){
  if(is.null(rewards)){
    stop('rewards should be a 1D numerical vector')
  }
  else {
    if(!is.numeric(rewards)) {
      if (class(rewards == 'matrix') & ncol != 1){
        stop('rewards should be a 1D numerical vector')
      }
      stop('rewards should be a 1D numerical vector')
    }
  }
  ######### Computation of T*, alpha and defect ##########
  rew_transformed = rewardtransformparm(rewards, obj$init_probs, obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(subint_mat = P, init_probs = alpha, defect = defect)
  attr(value, 'class') <- c('disc_phase_type')
  value
}

#' @export

segsites <- function(n, theta = 2){
  if (n<=1 | !is.numeric(n) | n %% 1 != 0 | n %% 1 != 0) {
    stop('n should be a positive integer larger than 1')
  }
  else if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }

  ph = cont_phase_type('T_Total', n = n)

  T_table = ph$subint_mat
  alpha = ph$init_probs
  defect = 0

  P = solve(diag(nrow(T_table)) - 2/theta * T_table)

  value <- list(subint_mat = T_table, init_probs = alpha, defect = 0)
  attr(value, "class") <- c("disc_phase_type")
  return(value)
}

#' @export

itons <- function(obj, itons, theta = 2) {
  if(!is.numeric(itons) | itons %% 1 != 0 | itons < 1){
    stop('itons should be a positive integer larger than 0')
  }
  else if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }

  if(length(obj$subint_mat) == 1){
    reward = matrix(obj$RewardM)
  }
  else{
    reward = obj$RewardM[,itons]
  }
  ######### Computation of T*, alpha and defect ##########
  rew_transformed = rewardtransformparm(reward, obj$init_probs, obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(subint_mat = P, init_probs = alpha, defect = defect)
  attr(value, 'class') <- c('disc_phase_type')
  value
}

#' @export

tailstat <- function(obj, k, theta = 2) {
  if(!is.numeric(k) | k %% 1 != 0 | k < 1){
    stop('k should be a positive integer larger than 0')
  }
  else if(!is.numeric(theta) | theta < 0){
    stop('theta should be a positive number')
  }

  if(nrow(obj$subint_mat) == k){
    reward = matrix(obj$RewardM[, ncol(obj$RewardM)])
  }
  else{
    reward = rowSums(obj$RewardM[, k:ncol(obj$RewardM)])
  }

  ######### Computation of T*, alpha and defect ##########
  rew_transformed = rewardtransformparm(reward, obj$init_probs, obj$subint_mat)
  alpha = rew_transformed$init_probs
  T_star = rew_transformed$subint_mat
  defect = rew_transformed$defect

  ########## Computation of P and p (transformation to DPH) ##########
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)

  value = list(subint_mat = P, init_probs = alpha, defect = defect)
  attr(value, 'class') <- c('disc_phase_type')
  value
}

