% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeMat.R
\name{TreeMat}
\alias{TreeMat}
\title{Create the tree matrix - Single Sample approach}
\usage{
TreeMat(orig.LT, orig.Y, orig.X)
}
\arguments{
\item{orig.LT}{Tree list}

\item{orig.Y}{original response vector}

\item{orig.X}{original covariate matrix}
}
\description{
This function generates the tree structure as a matrix, using the LTfunction for a single sample approach
}
\examples{
N<-2000
numx <- 5
alpha <-0.8
theta<-0.8
beta<- c(1,.8,.6,.4,.2)
gamma <- 1
Z <- rep(c(0, 1), each = N/2)
sigma <- diag(numx)
X.i <- mvrnorm(N,mu=rep(0,numx),Sigma=sigma)
W <- Z * ifelse(X.i[,1] > 0, 1, 0)
mu <- alpha + theta*Z + X.i \%*\% beta + W * gamma
Y <- rnorm(N, mean=mu)
matched2_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=X.i, version="prognostic")
LT <- LTfunction(Y, Z, X.i, X.i)
TreeMat(LT, matched2_orig$ymatch, as.matrix(matched2_orig$xmatch))
}
\keyword{matrix}
\keyword{tree}
