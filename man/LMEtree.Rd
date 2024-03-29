% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LMEtree.R
\name{LMEtree}
\alias{LMEtree}
\title{Create the tree as a list function}
\usage{
LMEtree(
  Y,
  X,
  ctl.ind,
  indx,
  indx.TEST,
  X.TEST,
  pval.thresh = 0.05,
  min.leaf.size = 12,
  met.type = "bonferroni",
  min.split.size = 7,
  split.pts = 0.05,
  sample = "single"
)
}
\arguments{
\item{Y}{response vector}

\item{X}{covariate matrix}

\item{ctl.ind}{treatment indicator}

\item{indx}{index}

\item{indx.TEST}{index in the test sample}

\item{X.TEST}{covariates correseponding to the test sample}

\item{pval.thresh}{p-value threshold}

\item{min.leaf.size}{minimum leaf size}

\item{met.type}{Default uses bonferroni correction. Other options include holm, hochberg, hommel, BH, BY, fdr and none. These are specified within the p.adjust function.}

\item{min.split.size}{minimum size to allow splitting}

\item{split.pts}{Split points}

\item{sample}{"single" or "double"}
}
\description{
This function uses the outcome (here the treatment effect) and covariates for the matched pairs and are determined from the MatchForTree function.
}
\keyword{Default}
\keyword{function}
\keyword{implementation}
\keyword{process}
\keyword{sample}
\keyword{single}
\keyword{splitting}
\keyword{the}
\keyword{tree}
\keyword{uses}
