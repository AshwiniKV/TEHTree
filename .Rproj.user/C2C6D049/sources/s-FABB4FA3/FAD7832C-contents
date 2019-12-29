#' Use the constructed single-sample tree and estimate the effects over the same sample
#'
#'@param orig.LT Use the constructed tree structure
#'@param orig.Y Use the response vector (original, not matched)
#'@param orig.X Use the covariate matrix (original, not matched)
#'@return A matrix composed of the tree structure and corresponding effects
#'
#'
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
