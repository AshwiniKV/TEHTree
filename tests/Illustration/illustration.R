
setwd("./Desktop/TEHTrees_code")
# Request access to the BLS dataset (6 months)
blsdata6<-read.csv("blsdata6.csv")
blsdata6<-blsdata6[,-1]

redbls<-blsdata6[,c(3, 10, 15, 53, 55, 84, 102, 106, 119:129, 132, 168, 172:174, 257:258, 273)]
redbls<-redbls[complete.cases(redbls),]

redbls<-redbls[!redbls$trt == "Control",]

Y<-redbls$kcal24h6
Z<-redbls$trt
Z<-factor(Z)
levels(Z)<-c("Control", "Trt", "Trt")
levels(Z)<-c("FALSE", "TRUE")
redbls$trt<-Z

indx <- sapply(redbls, is.factor)
redbls[indx] <- lapply(redbls[indx], function(x) as.numeric(x))

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

# TEHTree - Single sample approach
# Variations introduced by the matching process
#set.seed(15)
set.seed(60)
library(TEHTree)
X<-redbls2[,c(2, 4:5, 8)]
matched1_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=X, version="prognostic")
LT1<-LTfunction(Y, as.logical(Z), X, X)
TreeMat(LT1, matched1_orig$Y.match, as.matrix(matched1_orig$X.match))

##################################################

# Permuted
library(permute)
set.seed(250)
ids_X<-1:nrow(X)
newdata<-data.frame(cbind(ids_X, X))
nids_X<-shuffle(newdata$ids_X)
newX<-newdata[nids_X,-1]

matched2_orig <- MatchForTree(Y=Y, Z=as.logical(Z), X=newX, version="prognostic")
LT2<-LTfunction(Y, as.logical(Z), newX, newX)
TreeMat(LT2, matched2_orig$Y.match, as.matrix(matched2_orig$X.match))
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
set.seed(100)
cases <- sample(seq_len(nrow(redbls2)), round(nrow(redbls2) * .67))
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

plot_pres_htes <- function(conf.preds, conf.int = FALSE, z = 1.96) {
  if (is.null(conf.preds$predictions) || length(conf.preds$predictions) == 0)
    stop("conf.preds must include a column called 'predictions'")

  # Order the patients by the value of the treatment effect (Increasing)
  out <- ggplot(
    mapping = aes(
      x = rank(conf.preds$predictions),
      y = conf.preds$predictions
    )
  ) +geom_point() +  theme_light() + theme(axis.title.x = element_text( size=14, face="bold"),
                                           axis.title.y = element_text( size=14, face="bold"),
                                           axis.text.y = element_text( size=14, face="bold"),
                                           axis.text.x = element_text( size=14, face="bold")) +
    labs(x = "Patient (Test dataset)", y = "Estimated Treatment Effect - Ordered") +
     ylim(-200,600)

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
p1<-plot_pres_htes(preds, conf.int = TRUE)
p1
#ggsave("est_trteff_cf.pdf", p1)
