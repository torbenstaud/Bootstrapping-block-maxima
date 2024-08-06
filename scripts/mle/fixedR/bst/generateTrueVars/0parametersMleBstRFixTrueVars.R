#this script contains all the fixed parameters the other scripts use

alphaVec <- c(0.5, 1, 1.5) #frechet Shape in the limiting distribution
tsModVec <- c(3) # 1 =iid, 2 = m-dep, 3 = armax
distrVec <- c(3 ,4) #3=frechet, 4 = pareto
betaVec <- c(0, 0.25, 0.5, 0.75) # 0 = iid
mVec <- seq(40,100, by = 10) #numbers of effective sample sizes
#Fixed Parameters
beta <- 0.5 #armax dependence parameter
rVec <- c(90) #FIXED R
r <- rVec[1] #FIXED R
nVec <- rVec * mVec




