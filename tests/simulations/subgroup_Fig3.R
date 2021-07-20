#############################
# Simulation study
#############################

N<-100
#N<-2000 # For estimation related scenarios
numx <- 5
alpha<-0.8
theta<-0.8
beta<- c(1,.8,.6,.4,.2,0,0,0)
gamma <- c(0, 0, 0, 0) #might have to change this depending on scenario
phi <- c(0, 0,0, 0, 0) #might have to change this depending on scenario
rho<-0 #might have to change this depending on scenario
eta<-0
version="prognostic"
Y <- list()
X <- list()
trt.eff<-list()
X.i <- list()
W <- list()
M <- list()
mu <- list()

set.seed(205)

for(i in 1:1000){
  Z <- rep(c(0, 1), each = N/2)

  ### Define the correlation matrix
  sigma <- diag(numx)
  sigma[upper.tri(sigma) | lower.tri(sigma)] <- rho ## Compound symmetry

  X.i[[i]] <- mvrnorm(N,mu=rep(0,numx),Sigma=sigma)
  #X.i[[i]] <- matrix(rbinom(N*numx,size=1,prob=0.5), nrow=N)

  #might have to change this depending on simulation scenario
  X[[i]] <- cbind(X.i[[i]],
                  X.i[[i]][,1] * X.i[[i]][,2],
                  sqrt(abs(X.i[[i]][,3])),
                  as.integer(X.i[[i]][,4] > 0))

  #might have to change this depending on simulation scenario
  W[[i]] <- cbind(Z * ifelse(X[[i]][,1] > 0, 1, 0),
                  Z * X[[i]][,1],
                  Z * (sin(eta*X[[i]][,1])),
                  Z * (X[[i]][,1] > -0.5 & X[[i]][,1] < 0.5))

  #might have to change this depending on simulation scenario
  M[[i]] <- ifelse(X.i[[i]] > 0, 1, 0)

  mu[[i]] <- alpha + theta * Z + X[[i]] %*% beta + W[[i]] %*% gamma + M[[i]] %*% phi

  Y[[i]] <- rnorm(N, mean=mu[[i]])

  #might have to change this depending on simulation scenario
  trt.eff[[i]] <- theta + cbind(ifelse(X[[i]][,1] > 0, 1, 0),
                                ifelse(abs(X[[i]][,1]) < 0.8, 1, 0),
                                 ifelse(X[[i]][,1]>0 & X[[i]][,2]>0, 1, 0),
                                (X[[i]][,1] > -0.75 & X[[i]][,1] < 0.75)) %*% gamma
}

######################################################################################
# Causal Forest and Conditional inference tree combination
# Study the ability to extract groups

# Type I error
library(grf)
library(partykit)
cfobject<-mclapply(1:500, function(i) causal_forest(X.i[[i]], Y[[i]], Z), mc.cores = 5)

# Number of terminal nodes
extract<-function(cov, tree_obj){
  covariates<-data.frame(cov)
  effects<-as.numeric(tree_obj$predictions)
  newdata<-data.frame(effects, covariates)
  names(newdata)<-c("effects", "X1", "X2", "X3", "X4", "X5")
  trt_eff<-ctree(effects~X1+X2+X3+X4+X5, data = newdata)
  width(trt_eff)
}
trteff<-mclapply(1:500, function(i) extract(X.i[[i]], cfobject[[i]]), mc.cores = 5)
mean(unlist(trteff), na.rm = TRUE)

# X1 - Root node
extract<-function(cov, tree_obj){
  covariates<-data.frame(cov)
  effects<-as.numeric(tree_obj$predictions)
  newdata<-data.frame(effects, covariates)
  names(newdata)<-c("effects", "X1", "X2", "X3", "X4", "X5")
  trt_eff<-ctree(effects~X1+X2+X3+X4+X5, data = newdata)
  if(depth(trt_eff)>0){
  ifelse(strsplit(partykit:::.list.rules.party(trt_eff)[1], " ")[[1]][1] == "X1", 1, 0)
  }else{0}
  }
trteff<-mclapply(1:500, function(i) extract(X.i[[i]], cfobject[[i]]), mc.cores = 5)
mean(unlist(trteff), na.rm = TRUE)

# Any node- X1
extract<-function(cov, tree_obj){
  covariates<-data.frame(cov)
  effects<-as.numeric(tree_obj$predictions)
  newdata<-data.frame(effects, covariates)
  names(newdata)<-c("effects", "X1", "X2", "X3", "X4", "X5")
  trt_eff<-ctree(effects~X1+X2+X3+X4+X5, data = newdata)
  if(depth(trt_eff)>0){
  index<-strsplit(partykit:::.list.rules.party(trt_eff), " ")
  number<-function(i){
    ind<-seq(1, length(index[[i]]), 4)
    unlist(strsplit(index[[i]], " ")[ind])
  }
  vars<-unique(unlist(lapply(1:length(index), function(j) number(j))))
  ifelse("X1" %in% vars, 1, 0)
  }
  else{0}
}
trteff<-mclapply(1:500, function(i) extract(X.i[[i]], cfobject[[i]]), mc.cores = 5)
mean(unlist(trteff), na.rm = TRUE)


# Any split
extract<-function(cov, tree_obj){
  covariates<-data.frame(cov)
  effects<-as.numeric(tree_obj$predictions)
  newdata<-data.frame(effects, covariates)
  names(newdata)<-c("effects", "X1", "X2", "X3", "X4", "X5")
  trt_eff<-ctree(effects~X1+X2+X3+X4+X5, data = newdata)
  ifelse(depth(trt_eff)>0, 1, 0)
}

trteff<-mclapply(1:500, function(i) extract(X.i[[i]], cfobject[[i]]), mc.cores = 5)
mean(unlist(trteff), na.rm = TRUE)

#########################################################

# Causal Tree

library(causalTree)

# Terminal nodes
subgroup<-function(out, var, trt){
  newdata<-data.frame(out, var)
  names(newdata)<-c("Y", "X1", "X2", "X3", "X4", "X5")
  tree <- causalTree(out~., data = newdata, treatment = trt,
                     split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5,
                     cp = 0, minsize = 20, propensity = 0.5)
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  opfit <- prune(tree, opcp)
  nodes<-length(which(opfit$frame$var == "<leaf>"))
  nodes
}
subgrp<-mclapply(1:1000, function(i) subgroup(Y[[i]], X.i[[i]], Z), mc.cores = 5)
mean(unlist(subgrp))

# X1 - Root node
subgroup<-function(out, var, trt){
  newdata<-data.frame(out, var)
  names(newdata)<-c("Y", "X1", "X2", "X3", "X4", "X5")
  tree <- causalTree(out~., data = newdata, treatment = trt,
                     split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5,
                     cp = 0, minsize = 20, propensity = 0.5)
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  opfit <- prune(tree, opcp)
  #sub<-ifelse(opfit$frame$var[1] != "<leaf>", 1, 0)
  #sub<-ifelse(opfit$frame$var[1] == "X1"| opfit$frame$var[1] == "X2", 1, 0)
  sub<-ifelse(opfit$frame$var[1] == "X1", 1, 0)
  sub
}
subgrp<-mclapply(1:1000, function(i) subgroup(Y[[i]], X.i[[i]], Z), mc.cores = 5)
mean(unlist(subgrp))

# X1 - Any node

subgroup<-function(out, var, trt){
  newdata<-data.frame(out, var)
  names(newdata)<-c("Y", "X1", "X2", "X3", "X4", "X5")
  tree <- causalTree(out~., data = newdata, treatment = trt,
                     split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5,
                     cp = 0, minsize = 20, propensity = 0.5)
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  opfit <- prune(tree, opcp)
  sub<-ifelse(any(which(opfit$frame$var=="X1")), 1, 0)
  sub
}
subgrp<-mclapply(1:1000, function(i) subgroup(Y[[i]], X.i[[i]], Z), mc.cores = 5)
mean(unlist(subgrp))


# Any split
subgroup<-function(out, var, trt){
  newdata<-data.frame(out, var)
  names(newdata)<-c("Y", "X1", "X2", "X3", "X4", "X5")
  tree <- causalTree(out~., data = newdata, treatment = trt,
                     split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, xval = 5,
                     cp = 0, minsize = 20, propensity = 0.5)
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  opfit <- prune(tree, opcp)
  sub<-ifelse(opfit$frame$var[1] != "<leaf>", 1, 0)
  #sub<-ifelse(opfit$frame$var[1] == "X1"| opfit$frame$var[1] == "X2", 1, 0)
  #sub<-ifelse(opfit$frame$var[1] == "X1", 1, 0)
  sub
}
subgrp<-mclapply(1:1000, function(i) subgroup(Y[[i]], X.i[[i]], Z), mc.cores = 5)
mean(unlist(subgrp))

# TEHTree

# Estimation Error
# Double sample tree
Nused <- N/2
subjects <- c(1:(Nused/2), 151:(150+Nused/2))

#this gets the tree for the estimation sample
res <- lapply(1:500, function(i) LTfunction(Y[[i]][subjects], Z[subjects], X[[i]][subjects,], X.i[[i]][subjects,]), sample = "double")

predvar<-function(i){
  ifelse(!is.na(res[[i]]$splitleft),strsplit(res[[i]]$splitleft, " ")[[1]][1] == "X1", FALSE)}
preds_var<-lapply(1:1000, function(i) predvar(i))
mean(unlist(preds_var))

#######################################################
# Single Sample tree

ttree<-function(out, cov, trt){
  matched <- MatchForTree(Y=out,
                          Z = trt,
                          X=cov,
                          version="prognostic")
  ymatch <- matched$Y.match
  xmatch <- matched$X.match
  itrt <- matched$itrt
  ictl <- matched$ictl

  LT<-LMEtree(ymatch,
              as.matrix(xmatch),
              ictl,
              1:nrow(as.matrix(xmatch)), pval.thresh = 0.05, sample = "single")

  TreeMat(LT, ymatch, as.matrix(xmatch))
}

ttree_res<-mclapply(c(1:500), function(i) ttree(Y[[i]], X.i[[i]], Z), mc.cores = 5)

# Terminal nodes
predvar<-function(i){
  length(which(ttree_res[[i]]$split == "<leaf>"))
}
preds_var<-lapply(1:500, function(i) predvar(i))
mean(unlist(preds_var))

# Covariate X1 - Root node
predvar<-function(i){
  ifelse(!ttree_res[[i]]$split[1] == "<leaf>",
         strsplit(as.character(ttree_res[[i]]$split[1]), " ")[[1]][1] == "X1", 0)
}
preds_var<-lapply(1:500, function(i) predvar(i))
mean(unlist(preds_var))


# Covariate X1 - Any node
predvar<-function(i){
  if(!ttree_res[[i]]$split[1] == "<leaf>"){
index<-which(!ttree_res[[i]]$split == "<leaf>")
char<-as.character(ttree_res[[i]]$split[index])
numb<-seq(1, 3*length(index), 3)
term<-unlist(strsplit(char, " "))[numb]
}
ifelse(!ttree_res[[i]]$split[1] == "<leaf>", any(which(term == "X1")), 0)
}
preds_var<-lapply(1:500, function(i) predvar(i))
mean(unlist(preds_var))

# Presence of a split
predvar<-function(i){
  ifelse(!ttree_res[[i]]$split[1] == "<leaf>", 1, 0)
}
preds_var<-lapply(1:200, function(i) predvar(i))
mean(unlist(preds_var))

