#' Matching Function - Single and Double Sample Approach
#'
#' This function is used for matching with replacement and determines the treatment effect and relevant covariates for the matched pairs.
#'
#' @param Y response vector
#' @param X covariate matrix
#' @param Z treatment indicator
#' @param version version type of matching to be performed
#' @keywords matching prognostic treatment
#' @export
#' @import SuperLearner
#' @import Matching
#' @importFrom MASS mvrnorm
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
#' MatchForTree(Y, Z, X.i)

MatchForTree <- function(Y, Z, X, version = "prognostic") {

  if(is.null(dim(X))) {
    X.cat <- length(table(X)) <= 5
  }  else {
    X.cat <- apply(X,2,function(x) { length(table(x))<=5})
  }

  if(version == "standard") {
    M <- Match(Y=Y,
               Tr=Z,
               X=X,
               exact=X.cat,
               ties=FALSE,
               estimand="ATT",
               version = "fast")
    Yout <- Y[M$index.treated] - Y[M$index.control]

    if(is.null(dim(X))) {
      ##Xout <- (X[M$index.treated] + X[M$index.control]) / 2
      Xout <- min(abs(X[M$index.treated]), abs() + X[M$index.control]) / 2
    } else{
      Xout <- sapply(1:ncol(X),function(j) {
        x <- X[,j]
        if(X.cat[j]) { return(x[M$index.treated]) }
        else { return((x[M$index.treated]+x[M$index.control])/2)}
      })
    }

  }

  if(version == "prognostic") {
    #rhs <- paste( sprintf("X[, %d]", 1:(ncol(X))), collapse = " + ")
    #fmla <- sprintf("Y ~ %s", rhs)

    #D <- data.frame( Y, X )
    PS <- predict( SuperLearner( Y = Y, X = data.frame(X),
                                 SL.library = c("SL.glm", "SL.randomForest", "SL.gam",
                                                "SL.polymars", "SL.mean", "SL.glm.interaction",
                                                "SL.step", "SL.step.interaction")),
                   newdata = data.frame(X))$pred

    #    PS <- predict( lm( Y ~ .*., data = D ) )

    M <- Match(Y=Y,
               Tr=Z,
               X=PS,
               ties=FALSE,
               estimand="ATT",
               version = "fast")
    Yout <- Y[M$index.treated] - Y[M$index.control]

    if(is.null(dim(X))) {
      Xout <- X[M$index.treated]
      #Xout <- (X[M$index.treated] + X[M$index.control]) / 2
    } else{
      Xout <- X[M$index.treated, ]
      #Xout <- sapply(1:ncol(X),function(j) {
      #  x <- X[,j]
      #  if(X.cat[j]) { return(x[M$index.treated]) }
      #  else { return((x[M$index.treated]+x[M$index.control])/2)}
      #})
    }
  }

  if(version == "Xmin") {
    M <- Match(Y=Y,
               Tr=Z,
               X=X,
               exact=X.cat,
               ties=FALSE,
               estimand="ATT",
               version = "fast")
    Yout <- Y[M$index.treated] - Y[M$index.control]

    if(is.null(dim(X))) {
      ##Xout <- (X[M$index.treated] + X[M$index.control]) / 2
      XTC <- cbind( X[M$index.treated], X[M$index.control])
      Xout <- apply( XTC, 1, function(xrow) { xrow[which.min(abs(xrow))] })
    } else{
      Xout <- sapply(1:ncol(X),function(j) {
        x <- X[,j]
        if(X.cat[j]) { return(x[M$index.treated]) }
        else {
          XTC <- cbind( x[M$index.treated], x[M$index.control])
          return( apply( XTC, 1, function(xrow) { xrow[which.min(abs(xrow))] }) )
        } })  }
  }

  if(version == "Xtrt") {
    M <- Match(Y=Y,
               Tr=Z,
               X=X,
               exact=X.cat,
               ties=FALSE,
               estimand="ATT",
               version = "fast")
    Yout <- Y[M$index.treated] - Y[M$index.control]

    if(is.null(dim(X))) {
      ##Xout <- (X[M$index.treated] + X[M$index.control]) / 2
      Xout <- X[M$index.treated]
    } else{
      Xout <- sapply(1:ncol(X),function(j) {
        x <- X[,j]
        if(X.cat[j]) { return(x[M$index.treated]) }
        else {
          XTC <- cbind( x[M$index.treated], x[M$index.control])
          return( apply( XTC, 1, function(xrow) { xrow[which.min(abs(xrow))] }) )
        } })  }
  }

  if(version == "biasadjusted") {
    M <- Match(Y=Y,
               Tr=Z,
               X=X,
               exact=X.cat,
               ties=TRUE,
               estimand="ATT",
               BiasAdjust = TRUE,
               version = "standard")
    bias <- M$est.noadj - M$est
    Yout <- Y[M$index.treated] - Y[M$index.control] - bias

    if(is.null(dim(X))) {
      Xout <- (X[M$index.treated] + X[M$index.control]) / 2
    } else{
      Xout <- sapply(1:ncol(X),function(j) {
        x <- X[,j]
        if(X.cat[j]) { return(x[M$index.treated]) }
        else { return((x[M$index.treated]+x[M$index.control])/2)}
      })
    }

  }

  if(version == "autobalance") {
    sink("NUL")
    G <- GenMatch(Tr=Z, X=X, exact = X.cat, ties = FALSE, verbose = FALSE)
    sink()

    M <- Match(Y=Y,
               Tr=Z,
               X=X,
               exact=X.cat,
               ties=FALSE,
               estimand="ATT",
               Weight.matrix = G,
               version = "fast")
    Yout <- Y[M$index.treated] - Y[M$index.control]

    if(is.null(dim(X))) {
      Xout <- (X[M$index.treated] + X[M$index.control]) / 2
    } else{
      Xout <- sapply(1:ncol(X),function(j) {
        x <- X[,j]
        if(X.cat[j]) { return(x[M$index.treated]) }
        else { return((x[M$index.treated]+x[M$index.control])/2)}
      })
    }
  }

  return(list(Y.match=Yout,X.match=Xout, itrt = M$index.treated, ictl = M$index.control))
}
