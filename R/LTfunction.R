#' Determines the splitting process and split criterions
#'
#' This function generates the tree structure as a list function.
#'
#' @param Y response vector
#' @param Z treatment indicator
#' @param X covariate matrix
#' @param X.i covariate matrix - ids
#' @param version type of score - default is the prognostic score
#' @param pval.thresh Threshold of p-value
#' @param min.split.size Minimum size to split
#' @param sample Use a single or double sample - default is single. The double sample approach uses a separate sample to construct the tree and to estimate the effects. This can be modified as required.
#' @export
#' @importFrom stats rnorm as.formula p.adjust predict quantile sd
#' @importFrom MASS mvrnorm
#' @importFrom nlme lme
#' @examples
#' set.seed(10)
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
#' # Single sample approach
#' LT1<-LTfunction(Y, as.logical(Z), X.i, X.i)
#'
#' # Double sample approach
#' Nused <- 200
#' subjects <- c(1:(Nused/2), 1001:(1000+Nused/2))
#' LT2<-LTfunction(Y[subjects], Z[subjects], X.i[subjects,], X.i[subjects,], sample = "double")
#'
#'
LTfunction <- function(Y, Z, X, X.i, version = "prognostic", pval.thresh = 0.05,
                       min.split.size = 10, sample = "single"){

  if(sample == "single"){
    mat_orig <- MatchForTree(Y=Y, Z=as.logical(Z),
                             X=X.i, version="prognostic")
    matched<-mat_orig
    ymatch <- matched$Y.match
    xmatch <- matched$X.match
    itrt <- matched$itrt
    ictl <- matched$ictl
    LT <- LMEtree(ymatch, as.matrix(xmatch), ictl, 1:nrow(as.matrix(xmatch)), pval.thresh = 0.05,
                  min.split.size = min.split.size, sample = "single")
  }

  if(sample == "double"){
  n1 <- length(Y)/4
  n2 <- length(Y)/2
  n3 <- n1+n2

  matched <- MatchForTree(Y=Y[c(1:n1, (n2+1):n3)],
                          Z=Z[c(1:n1, (n2+1):n3)],
                          X=X[c(1:n1, (n2+1):n3),1:ncol(X.i)],
                          version=version)

  ymatch <- matched$Y.match
  xmatch <- matched$X.match
  itrt <- matched$itrt
  ictl <- matched$ictl
  #Y,X,ctl.ind,indx,indx.TEST,X.TEST,
  LT <- LMEtree(ymatch, as.matrix(xmatch), ictl,
                1:nrow(as.matrix(xmatch)),1:length(c((n1+1):n2,(n3+1):length(Y))),
                as.matrix(X[c((n1+1):n2,(n3+1):length(Y)),1:ncol(X.i)]), sample = "double", pval.thresh = 0.05)
  }

  return(LT)
}
