#'Build a single-sample tree
#'
#'@param Y Response vector
#'@param X Covariate matrix
#'@param ctl.ind Index corresponding to the control subject
#'@param indx Index to be selected for use in the tree-building process
#'@param pval.thresh Threshold of selected p-value
#'@param min.leaf.size Minium size of a terminal node
#'@param min.split.size Minimum size of a node to split
#'@param split.pts Set the increment in a sequence to search for split points
#'
#'@return A list object
#'
LMEstree <- function(Y,X,ctl.ind,indx,pval.thresh=0.2, min.leaf.size=12, min.split.size=7, split.pts=0.05,
                    min.eff.size=0.25*sd(Y,na.rm=T)) {

  leaf.list <- list(nodetype = "<leaf>",
                    index = indx,
                    splitleft = NA,
                    splitright = NA,
                    left = NA,
                    right = NA)

  #partysplit(length, index = indx)

  if(length(indx) < min.leaf.size) return(leaf.list)

  df  <- data.frame(Y,X,ctl.ind)[indx,]
  colnames(df) <- c("Y",paste0("X",1:ncol(X)),"ctl.ind")

  p <- sapply(1:ncol(X), function(j) {
    if(sum(diff(df[,sprintf("X%d",j)]))==0) { return(1) }
    fmla <- as.formula( sprintf("Y~X%d",j))
    fit  <- suppressWarnings(lme(fmla, random = ~1 | ctl.ind, data = df))
    ## Grab the p-value
    pval <- summary(fit)$tTable[2,5]
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
    split.points <- unique(quantile(df$X.split, seq(0,1,by=split.pts)))

    valid.split.points <- split.points[sapply(split.points, function(sp) {
      split.left <- indx[ which( df$X.split <= sp) ]
      split.right <- indx[ which( df$X.split > sp) ]
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
    ##if(abs(max(abs(coeffs))) < min.eff.size) return(leaf.list)
    if(!is.numeric(coeffs)) return(leaf.list)

    ## Find the largest coefficient
    wmax.coeff <- which.max(abs(coeffs))
    ## Split indx according
    split.left <- indx[ which( df$X.split < valid.split.points[wmax.coeff]) ]
    split.right <- indx[ which( df$X.split >= valid.split.points[wmax.coeff]) ]

    splitvarleft <- sprintf("X%d <= %.2f",k, valid.split.points[wmax.coeff])
    splitvarright <- sprintf("X%d > %.2f", k, valid.split.points[wmax.coeff])

    ## return split index
    return( list( nodetype = "internal",
                  index = indx,
                  splitleft = splitvarleft,
                  splitright = splitvarright,
                  left = LMEstree(Y,X,ctl.ind,split.left,pval.thresh),
                  right = LMEstree(Y,X,ctl.ind,split.right,pval.thresh) ) )
  }
  else {
    return( leaf.list ) ## STOP
  }

}
