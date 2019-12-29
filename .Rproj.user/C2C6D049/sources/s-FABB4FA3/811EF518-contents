#' Build the double-sample tree structure and specify the appropriate parameters
#'
#' @param Y Response vector
#' @param X Covariates matrix
#' @param ctl.ind Control Index
#' @param indx Index
#' @param indx.TEST Index in the Test dataset
#' @param X.TEST The defined covariate test dataset
#' @param pval.thresh p-value threshold
#' @param min.leaf.size Minimum lead size
#' @param min.split.size Minimum split size
#' @param split.pts Increment of the sequence to define split points
#' @export
#'
#'

LMEdtree <- function(Y,X,ctl.ind,indx,indx.TEST,X.TEST,pval.thresh=0.05,min.leaf.size=20,
                    min.split.size=10,split.pts=0.05) {

  leaf.list <- list(nodetype = "<leaf>",
                    index = indx,
                    index.TEST=indx.TEST,
                    splitleft = NA,
                    splitright = NA,
                    left = NA,
                    right = NA)

  if(length(indx) < min.leaf.size) return(leaf.list)

  df  <- data.frame(Y,X,ctl.ind)[indx,]

  df.TEST <- data.frame(Y,X.TEST,ctl.ind)[indx.TEST,]

  colnames(df) <- c("Y",paste0("X",1:ncol(X)),"ctl.ind")

  colnames(df.TEST) <- c("Y",paste0("X.TEST",1:ncol(X.TEST)),"ctl.ind")

  p <- sapply(1:ncol(X), function(j) {
    fmla <- as.formula( sprintf("Y~X%d",j))
    fit <- tryCatch(suppressWarnings(lme(fmla, random = ~1 | ctl.ind, data = df)),error=function(i){return(NA)})
    if(is.na(fit[1])){
      pval <- 1
    }else{
      ## Grab the p-value
      pval <- summary(fit)$tTable[2,5]
    }
    return(ifelse(is.nan(pval),1,pval))
  })

  if(is.null(p)) return(leaf.list)

  ## Apply the Bonferroni correction
  corr.pval <-  p.adjust(p, method = "bonferroni", n = length(p))

  k <- which.min(corr.pval)

  if(is.null(k)) return( leaf.list )
  if(corr.pval[k] < pval.thresh) {
    ## Do the splitting:
    df$X.split <- X[indx,k]

    df.TEST$X.split <- X.TEST[indx.TEST,k]

    split.points <- unique(quantile(df$X.split, seq(0,1,by=split.pts)))

    valid.split.points <- split.points[sapply(split.points, function(sp) {
      split.left <- indx[ which( df$X.split <= sp) ]
      split.right <- indx[ which( df$X.split > sp) ]

      split.left.TEST <- indx.TEST[ which( df.TEST$X.split <= sp) ]
      split.right.TEST <- indx.TEST[ which( df.TEST$X.split > sp) ]

      return( (length(split.left) >= min.split.size) & (length(split.right) >= min.split.size) )
    })]

    if(length(valid.split.points)==0) return(leaf.list) ## If no valid split points

    coeffs <- try( sapply( valid.split.points,
                           function(x) {
                             df$x <- x
                             fit <- try( lme(Y~I(X.split <= x), random = ~1 | ctl.ind, data = df) )
                             if(inherits(fit, "try-error")) return(0)
                             return(summary(fit)$tTable[2,1])
                           }) )
    if(inherits(coeffs, "try-error")) {
      warning("Error in splitting, cutting off tree building.")
      return(leaf.list)
    }

    ## Check size of coefficient
    if(!is.numeric(coeffs)) return(leaf.list)

    ## Find the largest coefficient
    wmax.coeff <- which.max(abs(coeffs))
    ## Split indx according
    split.left <- indx[ which( df$X.split <= valid.split.points[wmax.coeff]) ]
    split.right <- indx[ which( df$X.split > valid.split.points[wmax.coeff]) ]

    split.left.TEST <- indx.TEST[ which( df.TEST$X.split <= valid.split.points[wmax.coeff]) ]
    split.right.TEST <- indx.TEST[ which( df.TEST$X.split > valid.split.points[wmax.coeff]) ]

    splitvarleft <- sprintf("X%d <= %.2f",k, valid.split.points[wmax.coeff])
    splitvarright <- sprintf("X%d > %.2f", k, valid.split.points[wmax.coeff])

    splitvarleft.TEST <- sprintf("X.TEST%d <= %.2f",k, valid.split.points[wmax.coeff])
    splitvarright.TEST <- sprintf("X.TEST%d > %.2f", k, valid.split.points[wmax.coeff])

    ## return split index
    return( list( nodetype = "internal",
                  index = indx,
                  index.TEST = indx.TEST,
                  splitleft = splitvarleft,
                  splitright = splitvarright,
                  left = LMEtree(Y,X,ctl.ind,split.left, split.left.TEST,X.TEST,pval.thresh),
                  right = LMEtree(Y,X,ctl.ind,split.right, split.right.TEST,X.TEST,pval.thresh) ) )
  }
  else {
    return( leaf.list ) ## STOP
  }

}
