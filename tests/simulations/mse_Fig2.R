###################################One simulation#########################################
#Type I Error sims
N<-2000
numx <- 5
alpha<-0.8
theta<-0.8
beta<- c(1,.8,.6,.4,.2,0,0,0)
gamma <- c(0, 0, 0, 0) #might have to change this depending on scenario
phi <- c(0, 0,0, 0, 0) #might have to change this depending on scenario
rho<-0 #might have to change this depending on scenario
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
  W[[i]] <- cbind(Z * X[[i]][,1] * X[[i]][,2],
                  Z * ifelse(X[[i]][,2] > 0, 1, 0),
                  Z * (X[[i]][,1] > -0.5 & X[[i]][,1] < 0.5),
                  Z * ifelse(X[[i]][,2] > 0, 1, 0))

  #might have to change this depending on simulation scenario
  M[[i]] <- ifelse(X.i[[i]] > 0, 1, 0)

  mu[[i]] <- alpha + theta * Z + X[[i]] %*% beta + W[[i]] %*% gamma + M[[i]] %*% phi

  Y[[i]] <- rnorm(N, mean=mu[[i]])

  #might have to change this depending on simulation scenario
  trt.eff[[i]] <- theta + cbind(X[[i]][,1] * X[[i]][,2],
                                ifelse(X[[i]][,2] > 0, 1, 0),
                                (X[[i]][,1] > -0.5 & X[[i]][,1] < 0.5),
                                ifelse(X[[i]][,2] > 0, 1, 0)) %*% gamma


}

#might have to change this depending on simulation scenario
Nused <- 200
subjects <- c(1:(Nused/2), 1001:(1000+Nused/2))


#this gets the tree for the estimation sample
res <- lapply(1:1000, function(i) LTfunction(Y[[i]][subjects], Z[subjects], X[[i]][subjects,],
                                             X.i[[i]][subjects,], sample = "double"))

n1 <- (Nused/4)+1
n2<-(Nused/2)
n3 <- (1000+(Nused/4))+1
n4<-(1000+(Nused/2))

test.subj <- c(n1:n2,n3:n4)

#get the id of the covariate vector for which we want to find the treatment effect
test.case <- unlist(lapply(1:1000, function(i) which(X[[i]][test.subj,1] > 0 & Z[test.subj] == 1)[1]))

#get the ids in the leaf that the test case belongs to
grp.ids <- lapply(1:1000, function(i) grp_fnnc(res[[i]], test.case[i]))

#this gets the estimated treatment effect for the test case
#(and also for all covariate vectors that belong in the same leaf as the test case)
trt.effect <- lapply(1:1000, function(i) est_fncn(Y[[i]][test.subj], Z[test.subj],
                                  X[[i]][test.subj,], X.i[[i]][test.subj,], grp.ids[[i]], res[[i]]))
