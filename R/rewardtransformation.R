#' Functions used to transform to discrete phase type and calculate the number of segregating sites

#' rewardtransformation
#'
#' Function to compute reward transformed matrix T* as well as vector of initial probabilities
#' and the defect
#'
#' @param rewards reward vector
#' @param init_probs initial probabilities
#' @param subint_mat sub-intensity matrix
#'
#' @usage rewardtransformation(rewards, init_probs, subint_mat)
#'
#' @return A list containing `subint_mat` which is the reward transformed sub-intensity matrix, and `init_probs` which is the transformed initialization vector.
#'
#' @examples
#'
#' n = 8
#' obj = t_mrca(n)
#' reward_transformed_subint_matrix = rewardtransformation(c(3,2,1,0,0,1,0), obj$init_probs, obj$subint_mat)$subint_mat
#' print(reward_transformed_subint_matrix)
#'
#' @export

rewardtransformation = function(rewards, init_probs, subint_mat) {

  print("entered rewardtransformation")
  print(init_probs)
  print(rewards)
  print(subint_mat)


  # See BN p. 145

  # non-zero rewards
  d = sum(rewards > 0)

  # number of states
  p = length(init_probs)

  # transition matrix is -subint_ij minus subint_diagonal
  qmat = matrix(rep(0, p ^ 2), ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      #for (j in c(1:(i - 1), c((i + 1):p))) {
      if (i != j) {
        qmat[i, j] = -subint_mat[i, j] / subint_mat[i, i]
      }
    }
  }

  # If all rewards are positive
  if (d == p) {
    P_mat = qmat
    alpha_vec = init_probs
  }
  else{
    Qplusplus = qmat[(rewards > 0), (rewards > 0)]
    Qpluszero = qmat[(rewards > 0), (rewards == 0)]
    Qzeroplus = qmat[(rewards == 0), (rewards > 0)]
    Qzerozero = qmat[(rewards == 0), (rewards == 0)]

    # P_mat is the subintensity matrix of P_star_mat
    print("64")
    print(diag(1, nrow=p-d))
    print(Qzerozero)
    P_mat = Qplusplus + Qpluszero %*% solve(diag(1, nrow = p - d) - Qzerozero) %*% Qzeroplus

    piplus = init_probs[(rewards > 0)]
    pizero = init_probs[(rewards == 0)]

    # The initial distribution
    alpha_vec = piplus + pizero %*% solve(diag(1, nrow = p - d) - Qzerozero) %*% Qzeroplus

    subint_mat = as.matrix(subint_mat[(rewards > 0), (rewards > 0)])
    rewards = rewards[rewards > 0]
  }
  # Exit rate vector
  p_vec = 1 - rowSums(P_mat)
  T_star_mat = matrix(rep(0, d ^ 2), ncol = d)
  tstarvec = rep(0, d)

  for (i in 1:d) {
    # subintensity of Lambda_star_mat

    for (j in (1:d)[-i]) {
      T_star_mat[i, j] = -subint_mat[i, i] / rewards[i] * P_mat[i, j]
    }

    # .. and its exit rate vector
    tstarvec[i] = -subint_mat[i, i] / rewards[i] * p_vec[i]

    T_star_mat[i, i] = -sum(T_star_mat[i,]) - tstarvec[i]
  }

  cont_phase_type(subint_mat = T_star_mat,
                  init_probs = alpha_vec)
}
