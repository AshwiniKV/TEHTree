#' Estimation function over the test sample to determine the CATE (double-sample tree)
#'
#'@param Y Response vector
#'@param Z Treatment indicator
#'@param X Covariate dataset
#'@param X.i Covariate dataset corresponding to a test sample
#'@param grp.ids Determine the index values relevant to the identified subgroups
#'@param LT Tree structure
#'

est_fncn <- function(Y, Z, X, X.i, grp.ids, LT){

  n1 <- length(Y)/2
  ind1 <- grp.ids[which(grp.ids > n1)]
  ind2 <- grp.ids[which(grp.ids <= n1)]
  m1 <- mean(Y[ind1])
  m2 <- mean(Y[ind2])

  mean <- m1 - m2

  return(mean)
}
