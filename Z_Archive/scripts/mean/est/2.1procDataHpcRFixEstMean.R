########################################################
# Run file - Process Bootstrapped blockmaxima Frechet MLE for RFix
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


# gammaVec #frechet Shape in the limit
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



fullVarTibMeanEst <- tibble()
for(gammaInd in seq_along(gammaVec)){
  gamma <- gammaVec[gammaInd]
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      for(betaInd in seq_along(betaVec)){
        beta <- betaVec[betaInd]
        nameLoad <-paste(
          "MeanProcEst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
          "_beta", beta, ".Rdata", sep = ""
        )#name is: varTibMeanEst
        load(here(paste0("data/", nameLoad)))
        fullVarTibMeanEst <- bind_rows(
          fullVarTibMeanEst, varTibMeanEst
        )
      }
    }
  }
}

save(fullVarTibMeanEst,
     file = here(paste0("data/fullVarTibMeanEst")))