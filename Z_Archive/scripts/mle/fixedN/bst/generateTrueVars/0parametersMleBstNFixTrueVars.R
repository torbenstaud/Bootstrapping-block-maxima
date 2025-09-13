#this script contains all the fixed parameters the other scripts in this folder use

alphaVec <- c(0.5, 1, 1.5) #frechet Shape in the limiting distribution
tsModVec <- c(3) # 1 =iid, 2 = m-dep, 3 = armax
distrVec <- c(3 ,4) #3=frechet, 4 = pareto
betaVec <- c(0, 0.25, 0.5, 0.75) # 0 = iid
#Fixed Parameters
#N <- 5*10^3 #number of experiment repitions
#B <- 10^3 #number of Bootstrap replicates
nVec <- 10^3 #FIXED n
n <- nVec[1] #FIXED n
rVec <- seq(8, 40, by = 4)
mVec <- as.integer(n/rVec)
#kVec <- c(0, 1,2,3) #different k-max estimators
#dims <- c(N, length(mVec), length(kVec)) #shape of data array




