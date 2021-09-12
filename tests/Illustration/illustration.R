
# Request access to the BLS dataset (6 months)
library(TEHTree)
library(causalTree)
library(permute)
library(grf)

data("blsdata6")
names(blsdata6)<-c("trt", "bmi", "kcal24h6", "age", "edeq14", "hunger")
head(blsdata6)

# Causal Tree - double sample
# Use only 4 covariates - age, bmi, edeq14_0 and hunger
redbls1<-blsdata6[,-c(1)]
tree <- causalTree(kcal24h6~., data = redbls1, treatment = (blsdata6$trt -1),
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
X<-blsdata6[,c(2, 4:6)]
matched1_orig <- MatchForTree(Y=blsdata6$kcal24h6, Z=blsdata6$trt-1, X=X, version="prognostic")
LT1<-LTfunction(blsdata6$kcal24h6, blsdata6$trt-1, X, X)
TreeMat(LT1, matched1_orig$Y.match, as.matrix(matched1_orig$X.match))

##################################################

# Permuted - One example
# Variations introduced by the matching function and the permuted dataset

#set.seed(50)
ids_X<-1:nrow(X)
newdata<-data.frame(cbind(ids_X, X))
nids_X<-shuffle(newdata$ids_X)
newX<-newdata[nids_X,-1]

matched2_orig <- MatchForTree(Y=blsdata6$kcal24h6, Z=blsdata6$trt-1, X=newX, version="prognostic")
LT2<-LTfunction(blsdata6$kcal24h6, blsdata6$trt-1, newX, newX)
TreeMat(LT2, matched2_orig$Y.match, as.matrix(matched2_orig$X.match))
#####################################

newdata1<-data.frame(cbind(blsdata6$kcal24h6, newX))
names(newdata1)<-c("Y", "bmi0", "age", "edeq14.0", "hunger")
tree <- causalTree(Y~., data = newdata1, treatment = blsdata6$trt-1,
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F,
                   xval = 6, cp = 0, minsize = 10)

opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
rpart.plot(opfit)

###############
# Causal Forest
###############
set.seed(100)
cases <- sample(seq_len(nrow(blsdata6)), round(nrow(blsdata6) * .67))
train <- blsdata6[cases,]
test <- blsdata6[-cases,]

cf <- causal_forest(
  X = model.matrix(~ ., data = train[, c(2, 4:6)]),
  Y = train$kcal24h6,
  W = as.numeric(train$trt)-1,
  num.trees = 15000,
  seed = 10
)

preds <- predict(
  object = cf,
  newdata = model.matrix(~ ., data = test[, c(2, 4:6)]),
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
