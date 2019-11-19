#' Simulate Matrices from Kingsman coalescent
#'
#' Wrapper function that computes the matrices descrbining Kingsman coalescent ready for reward
#' transformation.
#'
#' @usage kingsman(n)
#'
#' @param n sample size (should be a positive integer)
#'
#' @return A \code{mult_phase_type} object containing the subintensity matrix, reward matrix and vector of initial probabilities
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
