library(partitions)

#' RewardTransformParm
#'
#' Function to compute reward transformed matrix T* as well as vector of initial probabilities and the defect
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
  list("newinitprob" = alphavec, "newsubintensitymatrix" = Tstarmat, "defect" = 1 - sum(alphavec))
}


#' RateMAndStateSpace
#'
#' Rate Matrix with corresponding state space
#'
#' @usage RateMAndStateSpace(n)

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
  return(list(RateM=RateM,StSpM=StSpM))
}


#' iton_mats
#'
#' Computes T* matrix, initial probabilities and defect for some frequency count (singleton, doubleton, ...)
#'
#' @usage iton_mats(n)
#'
#' for singletons
#' @usage iton_mats(n, itons = 1, theta = 3)

iton_mats <- function(n, init_probs = NA, itons = 0, theta = 2){

  ############## Step1: Preparation of Rate matrix (T), initial probabilities (pi) and reward vector ##############
  matrixes = RateMAndStateSpace(n)
  # Rate Matrix
  if(n == 2){
    T_table = matrix(matrixes$RateM[1, 1])
  }
  else{
    T_table = matrixes$RateM[1:ncol(matrixes$RateM)-1, 1:ncol(matrixes$RateM)-1]
  }

  # Initial Distribution (init_probs)
  # if nothing is supplied in the function argument, it creates a vector with first entry
  # being 1 and all the others 0

  if(is.na(init_probs)){
    if(n == 2){
      init_probs = matrix(c(1))
    }
    else{
      init_probs = c(1, rep(0, nrow(T_table)-1))
    }
  }

  # Specifying if we are considering all segregating sites or something more specific
  if(itons == 0){ # means all (singletons + doubletons + ...)
    if(n == 2){
      reward = matrix(sum(matrixes$StSpM[1,]))
    }
    else{
      reward = apply(matrixes$StSpM[1:nrow(matrixes$StSpM)-1, ], 1, sum)
    }
  }
  else{
    if(itons < n){
      if(n == 2){
        reward = matrix(matrixes$StSpM[1,1])
      }
      else{
        reward = matrixes$StSpM[1:nrow(matrixes$StSpM)-1,][,itons]
      }
    }
    else{
      return(0)
    }
  }
  ############## Step2: Computation of T*, alpha and defect ##############
  rew_transformed = rewardtransformparm(reward, init_probs, T_table)
  alpha = rew_transformed$newinitprob
  T_star = rew_transformed$newsubintensitymatrix

  list(T_star = T_star, alpha = alpha, defect = rew_transformed$defect)
}

# @describeIn dphtype
#'
#' pdphtype
#'
#' Distribution function.
#'
#' @usage pdphtype(k, alpha, T_star, theta)
#'
#' @export

dsegsites <- function(k, alpha, T_star, theta){
  P = solve(diag(nrow(T_star)) - 2/theta * T_star)
  p = matrix(1, nrow = nrow(P)) - P %*% matrix(1, nrow = nrow(P))

  dens_vec = rep(NA, k+1)
  for (i in 1:(k+1)){
    dens = alpha %*% (P%^%(i-1)) %*% p
    dens_vec[i] <- dens
  }
  dens_vec
}
