#' Create the tree matrix - Single Sample approach
#'
#' This function generates the tree structure as a matrix, using the LTfunction for a single sample approach
#'
#' @param orig.LT Tree list
#' @param orig.Y original response vector
#' @param orig.X original covariate matrix
#' @keywords tree matrix
#' @export
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
#' matched2_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=X.i, version="prognostic")
#' LT <- LTfunction(Y, Z, X.i, X.i)
#' TreeMat(LT, matched2_orig$ymatch, as.matrix(matched2_orig$xmatch))

TreeMat <- function(orig.LT, orig.Y, orig.X) {

  k <- 1
  frame.list <- list()

  TreeInfo <- function(LT, Y, X, k) {

    if(LT$nodetype == "<leaf>") {
      indx <- LT$index
      mn <- mean(Y[indx])
      frame.list[[k]] <<- c("<leaf>",
                            length(indx), # Number of observations
                            length(indx), # Weights of observations; assumed all 1
                            mean( (Y[indx] - mean(Y[indx]))^2), # deviance
                            mean(Y[indx]), #fitted value
                            0, 0, 0) # complexity, ncompete, nsurrogate
      return(k+1)
    } else {
      indx <- LT$index
      frame.list[[k]] <<- c( LT$splitleft,
                             length(indx),
                             length(indx),
                             mean( (Y[indx] - mean(Y[indx]))^2),
                             mean(Y[indx]),
                             0, 0, 0)
      k.new <- k+1
      frame.list[[k.new]] <<- c( LT$splitright,
                                 length(indx),
                                 length(indx),
                                 mean( (Y[indx] - mean(Y[indx]))^2),
                                 mean(Y[indx]),
                                 0, 0, 0)
      k.new <- k.new + 1
      k.new <- TreeInfo(LT$left, Y, X, k.new)
      k.new <- TreeInfo(LT$right, Y, X, k.new)
      return(k.new)
    }
  }

  TreeInfo(orig.LT, orig.Y, orig.X, 1)
  M <- do.call("rbind", frame.list)
  colnames(M) <- c("split", "nobs", "wtobs", "dev", "fitted", "complexity", "ncompete", "nsurrogate")
  return( data.frame(M) )
}
