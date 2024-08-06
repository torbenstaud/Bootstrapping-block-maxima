########################################################
# Run file - Process Bootstrapped blockmaxima mean statistics
########################################################
#!/usr/bin/env Rscript

#load libraries
library(tidyverse)
library(here)
library(ggplot2)
library(Rcpp)
#source code
Rcpp::sourceCpp("/gpfs/project/tosta102/bootstrap/0sourceC.cpp")  #c Functions
source("/gpfs/project/tosta102/bootstrap/0source.R") #R Functions
#load the parameters for the study----
source(here("0parameters.R"))
# gammaVec #frechet Shape in der Grenzverteilung
# tsModVec # 1 =iid, 2 = m-dep, 3 = armax
# distrVec  #1= GPD-gamma , 2= Normal, 3=frechet, 4 = pareto
# 
# N  #number of experiment repitions
# B  #number of bootstrap repetitions
# beta  #armax dependence parameter
# rVec  #value of r (fixed r setting!)
# mVec  #numbers of effective sample sizes
# kVec  #different k-max estimators

#Initialize HPC index from for loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#load different parameters and create grid for the HPC to loop over----
combs <- expand.grid(gammaVec, tsModVec, distrVec)
# initialize the parameters from the HPC loop ----
gamma <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
cat("Parameter gamma, ts_mod, distr= ", gamma, ts_mod, distr)

#Compare Bootstraps for the variance of blockmaxima----
niveau <- 0.05 #level for the confidence interval
#calculate the true variance of a block maximum
truth <- getTruthC(rVec[1], distr, ts_mod, gamma, 10^6, beta)
#variances for the respective block maxima estimators----
#initialize tibbles for the loop
varTibMeanEst <- tibble() #tibble for the estimated means

#tibble for the confidence intervals
cITibMeanEst <- tibble() 

#define the name of the data file to load
nameEst <- paste(
  "meanEst_", "gamma", gamma, "_Ts", ts_mod, 
  "_Distr", distr, "R.data", sep = ""
)
#load the array with simulated variances of block maxima
load(here(paste0("data/", nameEst))) #array is called estArrayMean
#name for saveFile: will save a file for each HPC array job
nameSave <- paste(
  "ProcMeanEst_", "gamma", gamma, "_Ts", ts_mod, 
  "_Distr", distr, "R.data", sep = ""
)

print(dim(estArrayMean))
#ratio of estimator variance vs mean over all bootstrap variances
ratios <- (estArrayMean[,,,1] %>% apply(c(2,3), var)) / 
  ((estArrayMean[,,,seq(2,B+1)] %>% 
      apply(c(1,2,3), var)) %>% apply(c(2,3), mean))
#bias, var, mse of respective estimators
mseFun <- function(xx){
  bias <- mean(xx) - truth
  var <- var(xx)
  return(c(bias, var, bias^2 + var))
}
estChar <- apply(estArrayMean[,,,1], c(2,3), mseFun)
tmpVarTibMeanEst <- tibble(
  gamma = gamma, tsModel = ts_mod, marginal = distr, 
  m = rep(mVec, length(kVec)), r = rVec[1], 
  k = as.factor(rep(kVec, each = length(mVec))),
  ratio = as.numeric(ratios), empBias = as.numeric(estChar[1,,]), 
  empVar = as.numeric(estChar[2,,]), empMSE = as.numeric(estChar[3,,])
)
varTibMeanEst <- bind_rows(varTibMeanEst, tmpVarTibMeanEst)




#confidence intervals
#checks for a fixed NSeed, whether true mean is covered by the CI
ciCont <- function(xx){
  lowBas <- 2*xx[1] - quantile(xx[seq(2,B+1)], 1-niveau/2)
  uppBas <- 2*xx[1] - quantile(xx[seq(2,B+1)], niveau/2)
  quantSym <- abs(
    quantile(xx[seq(2,B+1)], 1-niveau/2) - quantile(xx[seq(2,B+1)], niveau/2))/2
  lowSym <- xx[1] - quantSym
  uppSym <- xx[1] + quantSym
  return(
    unname(c(
      as.integer(lowBas <= truth && truth <= uppBas), #contained in basic
      as.integer(uppBas < truth), #underestimated in basic
      as.integer(lowBas > truth), #overestimated in basic
      uppBas -lowBas, #width of basic CI,
      as.integer(lowSym <= truth && truth <= uppSym), #contained in symmetric CI
      as.integer(uppSym < truth), #underestimated in symmetric CI
      as.integer(lowSym > truth), #overestimated in symmetric CI
      uppSym -lowSym #width of symmetric CI,
    ))
  )
}
#tibble with empirical coverage rates and CI width
tmpArr <- 
  apply(aperm(apply(estArrayMean, c(1,2,3), ciCont), c(2,3,4,1)), c(2,3,4), mean)

tmpCovTib <- tibble(
  gamma = gamma, tsModel = ts_mod, marginal = distr, 
  m = rep(mVec, 2*length(kVec)), r = rVec[1], 
  k = as.factor(rep(rep(kVec, each = length(mVec)),2)), 
  type = rep(c("basic", "symmetric"), each = length(kVec)*length(mVec)),
  coverage = c(as.numeric(tmpArr[,,1]),as.numeric(tmpArr[,,5])), 
  overEst = c(as.numeric(tmpArr[,,2]),as.numeric(tmpArr[,,6])),
  underEst = c(as.numeric(tmpArr[,,3]),as.numeric(tmpArr[,,7])),
  avgWidth = c(as.numeric(tmpArr[,,4]),as.numeric(tmpArr[,,8])),
)


save(tmpVarTibMeanEst, tmpCovTib,
     file = here(paste0("data/", nameSave)))