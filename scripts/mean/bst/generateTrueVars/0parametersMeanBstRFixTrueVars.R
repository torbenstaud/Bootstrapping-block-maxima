#this script contains all the fixed parameters the other scripts use

gammaVec <- c(-0.2, -0.1, 0, 0.1, 0.2) 
tsModVec <- c(3) # 1 =iid, 2 = m-dep, 3 = armax
distrVec <- c(1 ,5)
betaVec <- c(0, 0.25, 0.5, 0.75) # 0 = iid
mVec <- seq(40,100, by = 10) #numbers of effective sample sizes
#Fixed Parameters
rVec <- c(90) #FIXED R
r <- rVec[1] #FIXED R
nVec <- rVec * mVec




