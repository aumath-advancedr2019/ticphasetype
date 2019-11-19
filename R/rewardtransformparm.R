#' Functions used to transform to discrete phase type and calculate the number of segregating sites

#' RewardTransformParm
#'
#' Function to compute reward transformed matrix T* as well as vector of initial probabilities
#' and the defect
#'
#' @param rewards reward vector
#' @param init_probs initial probabilities
#' @param subint_mat sub-intensity matrix
#'
# @import dplyr
#'
#' @usage rewardtransformparm(rewards, init_probs, subint_mat)

rewardtransformparm <- function(rewards, init_probs, subint_mat){
  d <- sum(rewards > 0)
  p <- length(init_probs)
  qmat <- matrix(rep(0,p^2), ncol = p)
  for(i in 1:p){
    for(j in (1:p)[-i]){
      qmat[i,j] <- -subint_mat[i,j]/subint_mat[i,i]
    }
  }
  # If all rewards are stricly postive everything is simpler
  if(d == p){
    pmat <- qmat
    alphavec <- init_probs
  }
  else{
    qplusplus <- qmat[(rewards > 0),(rewards > 0)]
    qpluszero <- qmat[(rewards > 0),(rewards == 0)]
    qzeroplus <- qmat[(rewards == 0),(rewards > 0)]
    qzerozero <- qmat[(rewards == 0 ), (rewards == 0)]
    pmat <- qplusplus + qpluszero %*% solve(diag(1, nrow = p-d)-qzerozero) %*% qzeroplus
    piplus <- init_probs[(rewards > 0)]
    pizero <- init_probs[(rewards == 0)]
    alphavec <- piplus + pizero %*% solve(diag(1, nrow = p-d)-qzerozero) %*% qzeroplus
    subint_mat <- as.matrix(subint_mat[(rewards > 0), (rewards >0)])
    rewards <- rewards[rewards > 0]
  }
  pvec <- 1 - rowSums(pmat)
  Tstarmat <- matrix(rep(0,d^2), ncol = d)
  tstarvec <- rep(0,d)
  for(i in 1:d){
    for(j in (1:d)[-i]){
      Tstarmat[i,j] <- -subint_mat[i,i]/rewards[i]*pmat[i,j]
    }
    tstarvec[i] <- -subint_mat[i,i]/rewards[i]*pvec[i]
    Tstarmat[i,i] <- -sum(Tstarmat[i,])-tstarvec[i]
  }
  cont_phase_type(subint_mat = Tstarmat, init_probs = alphavec)
}
