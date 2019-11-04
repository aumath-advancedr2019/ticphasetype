#' Functions used to transform to discrete phase type and calculate the number of segregating sites

#' RewardTransformParm
#'
#' Function to compute reward transformed matrix T* as well as vector of initial probabilities
#' and the defect
#'
#' @import partitions
# @import dplyr
#'
#' @usage rewardtransformparm(n, init_probs, subint_mat)
#'
#' @export

rewardtransformparm <- function(rewards, init_probs, subint_mat){
  d <- sum(rewards > 0)
  p <- length(init_probs)
  qmat <- matrix(rep(0,p^2), ncol = p)
  for(i in 1:p){
    for(j in (1:p)[-i]){
      qmat[i,j] <- -subint_mat[i,j]/subint_mat[i,i]
    }
  }
  ##
  # If all rewards are stricly postive everything is simpler
  ##
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
  phase_type(subint_mat = Tstarmat, init_probs = alphavec)
}

#----------------------------------------------------------------------------------------------

#' RateMAndStateSpace
#'
#' Rate Matrix with corresponding state space
#'
#' @usage RateMAndStateSpace(n)
#'
#' @export

RateMAndStateSpace <- function(n){
  ##----------------------------------------------------
  ## Possible states
  ##----------------------------------------------------
  ## Size of the state space (number of states)
  nSt <- P(n)
  ## Definition of the state space
  StSpM <- matrix(ncol=n,nrow=nSt)
  ## Set of partitions of [n]
  x <- parts(n)
  ## Rewriting the partitions as (a1,...,an)
  for (i in 1:nSt) {
    st <- x[,i]
    StSpM[i,] <- tabulate(x[,i],nbins=n)
  }
  ## Reordering
  StSpM <- StSpM[order(rowSums(StSpM),decreasing=TRUE),]
  ## Because of this ordering we can't 'go back', i.e.
  ## below the diagonal the entries are always zero
  ##----------------------------------------------------
  ## Intensity matrix
  ##----------------------------------------------------
  RateM <- matrix(0,ncol=nSt,nrow=nSt)
  ## Algorithm for finding rates between states
  for (i in 1:(nSt-1)){
    for (j in (i+1):nSt){
      cvec <- StSpM[i,]-StSpM[j,]
      check1 <- sum(cvec[cvec>0])==2
      check2 <- sum(cvec[cvec<0])==-1
      if (check1 & check2){
        ## Size(s) of the block(s) and the corresponding rates
        tmp <- StSpM[i,which(cvec>0)]
        RateM[i,j] <- ifelse(length(tmp)==1,tmp*(tmp-1)/2,prod(tmp))
      }
    }
  }
  ## Diagonal part of the rate matrix
  for (i in 1:nSt){
    RateM[i,i] <- -sum(RateM[i,])
  }
  list(RateM=RateM,StSpM=StSpM)
}
