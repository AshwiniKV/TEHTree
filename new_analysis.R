getwd()
setwd("./Desktop")
# Request access to the BLS dataset (6 months)
blsdata6<-read.csv("blsdata6.csv")
blsdata6<-blsdata6[,-1]

redbls<-blsdata6[,c(3, 10, 15, 53, 55, 84, 102, 106, 119:129, 132, 168, 172:174, 257:258, 273)]
redbls<-redbls[complete.cases(redbls),]

# Remove subjects in the control group
redbls<-redbls[!redbls$trt == "Control",]

Y<-redbls$kcal24h6
Z<-redbls$trt
Z<-factor(Z)
levels(Z)<-c("Control", "Trt", "Trt")
levels(Z)<-c("FALSE", "TRUE")
redbls$trt<-Z

indx <- sapply(redbls, is.factor)
redbls[indx] <- lapply(redbls[indx], function(x) as.numeric(x))

names(redbls)
kcaltrt<-redbls$kcal24h6[Z == "TRUE"]
kcalctrl<-redbls$kcal24h6[Z == "FALSE"]
mean(kcaltrt) - mean(kcalctrl)

###############
# Full dataset - 25 covariates

# Double sample Causal Tree
library(causalTree)
redbls1<-redbls[,-1]
tree <- causalTree(kcal24h6~., data = redbls1, treatment = as.logical(Z),
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F,
                   xval = 6, cp = 0, minsize = 10)

opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
rpart.plot(opfit)

#######################
set.seed(10)
library(SuperLearner)
library(mvtnorm)
library(nlme)
library(Matching)
library(parallel)
library(partykit)
library(rgenoud)
X<-redbls[,c(2:5, 7:27)]
# Use single sample TEHTree code
matched1_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=X, version="prognostic")
#matched2<-matched1_orig
matched<-matched1_orig
ymatch <- matched$Y.match
xmatch <- matched$X.match
itrt <- matched$itrt
ictl <- matched$ictl
LT1 <- LMEstree(ymatch,
              as.matrix(xmatch),
              ictl,
              1:nrow(as.matrix(xmatch)), pval.thresh = 0.85, min.split.size = 5)

TreeMat(LT1, ymatch, as.matrix(xmatch))
##################################################

# Permuted
library(permute)
ids_X<-1:nrow(X)
newdata<-data.frame(cbind(ids_X, X))
nids_X<-shuffle(newdata$ids_X)
newX<-newdata[nids_X,-1]

# TEHTree
matched2_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=newX, version="prognostic")
matched<-matched2_orig
ymatch <- matched$Y.match
xmatch <- matched$X.match

itrt <- matched$itrt
ictl <- matched$ictl
LT <- LMEtree(ymatch,
              as.matrix(xmatch),
              ictl,
              1:nrow(as.matrix(xmatch)), pval.thresh = 0.8, min.split.size = 5)

TreeMat(LT, ymatch, as.matrix(xmatch))
#####################################

# Causal Tree
library(causalTree)

newdata1<-data.frame(cbind(Y, newX))

tree <- causalTree(Y~., data = newdata1, treatment = as.logical(Z),
                   split.Rule = "CT", cv.option = "CT", split.Honest = T,
                   cv.Honest = T, split.Bucket = F, minsize = 10)

opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
rpart.plot(opfit)

###############
# Causal Forest
###############
# Reproducible result note: Differences introduced by the training and test samples.

cases <- sample(seq_len(nrow(redbls)), round(nrow(redbls) * .75))
train <- redbls[cases,]
test <- redbls[-cases,]

cf <- causal_forest(
  X = model.matrix(~ ., data = train[, c(2:5, 7:ncol(train))]),
  Y = train$kcal24h6,
  W = as.numeric(train$trt)-1,
  num.trees = 15000,
  seed = 10
)

preds <- predict(
  object = cf,
  newdata = model.matrix(~ ., data = test[, c(2:5, 7:ncol(test))]),
  estimate.variance = TRUE
)

plot_pres_htes <- function(conf.preds, conf.int = FALSE, z = 1.96) {
  if (is.null(conf.preds$predictions) || length(conf.preds$predictions) == 0)
    stop("conf.preds must include a column called 'predictions'")

  # Order the patients by the value of the treatment effect (Increasing)
  out <- ggplot(
    mapping = aes(
      x = rank(conf.preds$predictions),
      y = conf.preds$predictions
    )
  ) + geom_point() +
    labs(x = "Patient (Test dataset)", y = "Estimated Treatment Effect - Ordered") +
    theme_light() +ylim(-200, 600)

  if (conf.int && length(conf.preds$variance.estimates) > 0) {
    out <- out +
      geom_errorbar(
        # Compute the lower and upper interval value
        mapping = aes(
          ymin = conf.preds$predictions + z * sqrt(conf.preds$variance.estimates),
          ymax = conf.preds$predictions - z * sqrt(conf.preds$variance.estimates)
        )
      )
  }
  # ggplot
  return(out)
}

library(ggplot2)
# Plot the effects with the confidence interval

plot_pres_htes(preds, conf.int = TRUE)

################################################################################################################
################################################################################################################

# Subset of covariates

redbls2<-blsdata6[,c(3, 15, 84, 168, 122, 257:258, 174, 273)]
redbls2<-redbls[,c(1, 3, 6, 21, 12, 25:26, 24, 27)]
###############
# Causal Tree - double sample
# Use only 4 covariates - age, bmi, hunger and edeq14.0

library(causalTree)
redbls1<-redbls2[,-1]
redbls1<-redbls1[,-c(5:6, 8)]
tree <- causalTree(kcal24h6~., data = redbls1, treatment = as.logical(Z),
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F,
                   xval = 6, cp = 0, minsize = 10)

opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
rpart.plot(opfit)

##############################################

# TEHTree
X<-redbls2[,c(2, 4:5, 8)]
matched1_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=X, version="prognostic")
#matched2<-matched1_orig
matched<-matched1_orig
ymatch <- matched$Y.match
xmatch <- matched$X.match
itrt <- matched$itrt
ictl <- matched$ictl
LT1 <- LMEstree(ymatch,
               as.matrix(xmatch),
               ictl,
               1:nrow(as.matrix(xmatch)), pval.thresh = 0.05, min.split.size = 10)

TreeMat(LT1, ymatch, as.matrix(xmatch))
##################################################

# Permuted
library(permute)
ids_X<-1:nrow(X)
newdata<-data.frame(cbind(ids_X, X))
nids_X<-shuffle(newdata$ids_X)
newX<-newdata[nids_X,-1]

matched2_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=newX, version="prognostic")
matched<-matched2_orig
ymatch <- matched$Y.match
xmatch <- matched$X.match

itrt <- matched$itrt
ictl <- matched$ictl
LT <- LMEstree(ymatch,
              as.matrix(xmatch),
              ictl,
              1:nrow(as.matrix(xmatch)), pval.thresh = 0.05, min.split.size = 10)

TreeMat(LT, ymatch, as.matrix(xmatch))
#####################################

library(causalTree)

newdata1<-data.frame(cbind(Y, newX))
names(newdata1)<-c("Y", "bmi0", "age", "edeq14.0", "hunger")
tree <- causalTree(Y~., data = newdata1, treatment = as.logical(Z),
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F,
                   xval = 6, cp = 0, minsize = 10)

opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
rpart.plot(opfit)

###############
# Causal Forest
###############

cases <- sample(seq_len(nrow(redbls2)), round(nrow(redbls2) * .75))
train <- redbls2[cases,]
test <- redbls2[-cases,]

library(grf)

cf <- causal_forest(
  X = model.matrix(~ ., data = train[, c(2, 4:5, 8)]),
  Y = train$kcal24h6,
  W = as.numeric(train$trt)-1,
  num.trees = 15000,
  seed = 10
)

preds <- predict(
  object = cf,
  newdata = model.matrix(~ ., data = test[, c(2, 4:5, 8)]),
  estimate.variance = TRUE
)


