#' Functions used to transform to discrete phase type and calculate the number of segregating sites

#' RewardTransformParm
#'
#' Function to compute reward transformed matrix T* as well as vector of initial probabilities and the defect
#'
#' @import partitions
#' @import dplyr
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

#----------------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------------

#'
#'
#' Computes P and T* matrix, initial probabilities and defect for some frequency count (singleton, doubleton, ...)
#'
#' @usage iton_mats(n)
#'
#' for singletons
#' @usage iton_mats(n, itons = 1, theta = 3)

discph = function(n, k = 0, pi_vec = NA, itons = 0, theta = 2, moment = F, tail_stat = F){

  # All Segregating Sites - This computes the P matrix faster than setting rewards to sum(ri)
  if(itons == 0){
    ph = phase_type('T_Total', n = n)

    T_table = ph$subint_mat
    alpha = ph$init_probs
    defect = 0

    P = solve(diag(nrow(T_table)) - 2/theta * T_table)
  }

  # Special Case (eg. singletons or tail statistic)
  else{
    # There is no reason to proceed if this condition is True
    if(itons >= n){
      if(moment){
        return(list(Moment_1 = 0, Moment_2 = 0, Variance = 0, defect = 0))
      }
      else{return(0)}
    }

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
    if(is.na(pi_vec)){
      if(n == 2){
        pi_vec = matrix(c(1))
      }
      else{
        pi_vec = c(1, rep(0, nrow(T_table)-1))
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
    rew_transformed = rewardtransformparm(reward, pi_vec, T_table)
    alpha = rew_transformed$newinitprob
    T_star = rew_transformed$newsubintensitymatrix
    defect = rew_transformed$defect

    ########## Step3: Computation of P and p (transformation to DPH) ##########
    P = solve(diag(nrow(T_star)) - 2/theta * T_star)
  }
  # Calculating The Moments
  if(moment | tail_stat){
    m1 = alpha %*% solve(diag(nrow(P)) - P) %*% matrix(1, nrow = nrow(P)) # Corollary 1.2.64
    m2 = 2 * alpha %*% P %*% (solve(diag(nrow(P)) - P) %*% solve(diag(nrow(P)) - P)) %*% matrix(1, nrow = nrow(P))# Theorem 1.2.69
    es = m1 - 1 + defect
    mvar = m2 - es^2 - es
    return(list(Moment_1 = m1[1], Moment_2 = m2[1], Variance = mvar[1], defect = defect))
  }
  # Calculating the Probability P(S = k)
  else{
    if(itons == 0){
      list(P = P, alpha = alpha, defect = defect)
    }
    else{
      list(P = P, T_star = T_star, alpha = alpha, defect = defect)
    }
  }
}


sfs <- function(n_vec, theta = 2){
  sfs_df <- tibble()
  # iterate through sample sizes given in n_vec
  for (i in 1:length(n_vec)){

    # calculate E[ksi_i] and Var[ksi_i] for every i < n
    for (iton in 1:(n_vec[i]-1)){
      m = discph(n_vec[i], itons = iton, moment = T, theta = theta)
      varksi = m$Variance
      ksi = m$Moment_1 - 1 + m$defect
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

dsegsites <- function(x, alpha, P){
  p = matrix(1, nrow = nrow(P)) - P %*% matrix(1, nrow = nrow(P))
  dens_vec = c()
  for (i in x){
    dens_vec <- c(dens_vec, alpha %*% (P%^%(i)) %*% p)
  }
  dens_vec
}

