#this script contains all the fixed parameters the other scripts use

gammaVec <- c(-0.2, -0.1, 0, 0.1, 0.2) #frechet Shape in der Grenzverteilung
tsModVec <- c(3) # 1 =iid, 2 = m-dep, 3 = armax
distrVec <- c(1 ,5) #3=frechet, 4 = pareto
betaVec <- c(0, 0.25, 0.5, 0.75) # 0 = iid
mVec <- seq(40,100, by = 10) #numbers of effective sample sizes
N <- 5*10^3 #number of experiment repitions
B <- 10^3  #number of bootstrap repetitions
# n <- 10^3 NOT NEEDED AS WE ARE IN THE FIXED r SCENARIO

rVec <- c(90) #value of r (fixed r setting!)
mVec <- seq(40,100, by = 10)  #numbers of effective sample sizes: for true: 50 to 125
kVec <- c(0, 1,2,3) #different k-max estimators
#dims <- c(N, length(mVec), length(kVec), B+1) #shape of data array