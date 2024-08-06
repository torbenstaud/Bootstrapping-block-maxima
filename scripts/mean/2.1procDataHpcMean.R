########################################################
# Run file - Combine processed Bootstrapped blockmaxima mean statistics
########################################################
#!/usr/bin/env Rscript

#load libraries
library(tidyverse)
library(here)
library(gridExtra)

#sourcing of relevant files----
# the following files need to be sourced: 
# -0parameters.R in the same folder as this R.file
# gammaVec #frechet Shape in der Grenzverteilung
# tsModVec # 1 =iid, 2 = m-dep, 3 = armax
# distrVec  #1= GPD-gamma , 2= Normal, 3=frechet, 4 = pareto, 5 = GEV-gamma
# 
# N  #number of experiment repitions
# B  #number of bootstrap repetitions
# beta  #armax dependence parameter
# rVec  #value of r (fixed r setting!)
# mVec  #numbers of effective sample sizes
# kVec  #different k-max estimators


varTibMeanEst <- tibble() #tibble for the variances
if(FALSE){
  
#initialize Arrays for the histogramm slices
estArraySlices <- 
  array(dim = c(length(gammaVec), length(tsModVec), length(distrVec),
                N, length(mVec), length(kVec)
  )) 
bstArraySlices <- 
  array(dim = c(length(gammaVec), length(tsModVec), length(distrVec),
                B, length(mVec), length(kVec)
  ))
}
#tibble for the confidence intervals
cITibMeanEst <- tibble() 


for(gammaInd in seq_along(gammaVec)){
  gamma <- gammaVec[gammaInd]
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      #define the name of the data file to load
      nameLoad <- paste(
        "ProcMeanEst_", "gamma", gamma, "_Ts", ts_mod, 
        "_Distr", distr, "R.data", sep = ""
      )
      #load the array with simulated variances of block maxima
      load(here(paste0("data/", nameLoad))) 
      #has 2 (three) variables inside: tmpVarTibMeanEst, tmpCovTib
      #update the varTib
      varTibMeanEst <- bind_rows(varTibMeanEst, tmpVarTibMeanEst)
      
      #update the CI tibble
      cITibMeanEst <- bind_rows(cITibMeanEst, tmpCovTib)
      
    }
  }
}
save(cITibMeanEst, varTibMeanEst,
     file = here("data/processedData.RData"))
