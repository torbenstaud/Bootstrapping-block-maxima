########################################################
# Run file - Process Estimated blockmaxima mean for RFix
########################################################
#!/usr/bin/env Rscript

#load libraries
library(tidyverse)
library(here)
library(ggplot2)
library(Rcpp)

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parametersRFix.R in the same folder as this R.file


# alphaVec #frechet Shape in the limit
# tsModVec # 1 =iid, 2 = m-dep, 3 = armax
# distrVec  #1= GPD-gamma , 2= Normal, 3=frechet, 4 = pareto
# betaVec #armax dependence paramaters
#
# N  #number of experiment repitions
# B  #number of bootstrap repetitions
# beta  #armax dependence parameter
# rVec  #value of r (fixed r setting!)
# mVec  #numbers of effective sample sizes
# kVec  #different k-max estimators
cat("Parameter r")

#Initialize HPC index from for loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#load different parameters and create grid for the HPC to loop over----
combs <- expand.grid(gammaVec, tsModVec, distrVec, betaVec)
# initialize the parameters from the HPC loop ----
gamma <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
#Compare estimation of frechet shape----

#mean for the respective block maxima estimators----
#initialize tibbles for the loop
varTibMeanEst <- tibble() #tibble for the estimated gamma


trueMean <- getTruthC(r, distr, ts_mod, gamma, 10^6, beta)
nameLoad <-paste(
  "MleEst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
  "_beta", beta, ".Rdata", sep = ""
)
load(here(paste0("data/", nameLoad)))
vars <- as.numeric(apply(meanArray, c(2,3), var))
biases <- as.numeric(apply(meanArray, c(2,3), mean)) - trueMean
updTib <-
  tibble(
    gamma = gamma, tsMod = ts_mod, marginal = distr, beta = beta,
    m = rep(mVec, length(kVec)), r = rVec[1], 
    k = as.factor(rep(kVec, each = length(mVec))),
    varEst = vars, biasEst = biases
  ) %>% mutate(mse = varEst + biasEst^2)
varTibMeanEst <- bind_rows(varTibMeanEst, updTib)

nameSave <- paste(
  "MeanProcEst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
  "_beta", beta, ".Rdata", sep = ""
)
save(varTibMeanEst,
     file = here(paste0("data/", nameSave)))
cat(paste0("\n We had r = ", r))