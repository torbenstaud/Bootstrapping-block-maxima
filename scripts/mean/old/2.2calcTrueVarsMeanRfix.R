########################################################
# Run file - Generate Bootstrapped mean of block maxima statistics
########################################################
#!/usr/bin/env Rscript
#necessary libraries
library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths 
library(tidyverse) #tibbles

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parameters.R in the same folder as this R.file

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
combs <- expand.grid(gammaVec, tsModVec, distrVec, betaVec, mVec)
# initialize the parameters from the HPC loop ----
gamma <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
m <- combs[index, 5]
cat("Parameter gamma, ts_mod, distr, beta, m = ", gamma, ts_mod, distr, beta, m)
r <- rVec[1] #fixed r for this study
n <- m*r
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________





trueTmpVarsTib <- 
  tibble(
    gamma = gamma, tsModel = ts_mod, marginal = distr,
    m = m, r = r, beta = beta, k = 0, 
    trueSlidEstVar = trueMeanEstVarC(m, r, distr, ts_mod, gamma, 10^6, beta)
  )

cat(paste("\n We had r = ", r, " and gamma, tsMod, distr, m = ",gamma,
           ts_mod, distr, m))
nameSave <- paste0("trueTmpVarsTib",gamma,ts_mod,distr, beta, m,r)
save(trueTmpVarsTib, file = here(paste0("data/r", r, "/", nameSave)))