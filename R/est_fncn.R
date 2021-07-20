#' Estimation function - Double Sample Approach
#'
#' This function estimates the treatment effect using a known tree structure
#'
#' @param Y response vector
#' @param Z treatment indicator
#' @param X covariate matrix
#' @param X.i covariate matrix
#' @param grp.ids IDs of groups
#' @param LT Tree list - double sample
#' @export
#' @keywords sample tree
#' @examples
#' N<-2000
#' numx <- 5
#' alpha <-0.8
#' theta<-0.8
#' beta<- c(1,.8,.6,.4,.2)
#' gamma <- 1
#' Z <- rep(c(0, 1), each = N/2)
#' sigma <- diag(numx)
#' X.i <- mvrnorm(N,mu=rep(0,numx),Sigma=sigma)
#' W <- Z * ifelse(X.i[,1] > 0, 1, 0)
#' mu <- alpha + theta*Z + X.i %*% beta + W * gamma
#' Y <- rnorm(N, mean=mu)
#'
#' Nused <- 200
#' subjects <- c(1:(Nused/2), 1001:(1000+Nused/2))
#' #Using a double sample approach
#' res <- LTfunction(Y[subjects], Z[subjects], X.i[subjects,], X.i[subjects,], sample= "double")
#' n1<-(Nused/4)+1
#' n2<-(Nused/2)
#' n3<-(1000+(Nused/4))+1
#' n4<-(1000+(Nused/2))
#' test.subj <- c(n1:n2,n3:n4)
#' # Case when a treated patient has characteristics defined by X1 > 0
#' test.case <- which(X.i[test.subj,1] > 0 & Z[test.subj] == 1)[1]
#' # Get the ids in the leaf that the test case belongs to
#' group.ids <- grp_func(res, test.case)
#' trt.effect <- est_fncn(Y[test.subj], Z[test.subj], X.i[test.subj,], X.i[test.subj,], group.ids, res)

est_fncn <- function(Y, Z, X, X.i, grp.ids, LT){

  n1 <- length(Y)/2
  ind1 <- grp.ids[which(grp.ids > n1)]
  ind2 <- grp.ids[which(grp.ids <= n1)]
  m1 <- mean(Y[ind1])
  m2 <- mean(Y[ind2])

  mean <- m1 - m2

  return(mean)
}
