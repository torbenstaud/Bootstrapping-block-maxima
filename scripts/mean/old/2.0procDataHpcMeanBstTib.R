########################################################
# Run file - Calculate Bootstrap Variances for the mean estimator
########################################################
#!/usr/bin/env Rscript

#load libraries
library(tidyverse)
library(here)
library(gridExtra)

#sourcing of relevant files----
# the following files need to be sourced: 
# -0parameters.R in the same folder as this R.file
mVec <- seq(50,100, by = 5)
# tsModVec # 1 =iid, 2 = m-dep, 3 = armax
# distrVec  #1= GPD-gamma , 2= Normal, 3=frechet, 4 = pareto, 5 = GEV-gamma
# 
# N  #number of experiment repitions
# B  #number of bootstrap repetitions
# beta  #armax dependence parameter
# rVec  #value of r (fixed r setting!)
# mVec  #numbers of effective sample sizes
# kVec  #different k-max estimators


bstTib <- tibble()

for(gammaInd in seq_along(gammaVec)){
  gamma <- gammaVec[gammaInd]
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      #define the name of the data file to load
      nameLoad <- paste(
        "meanEst_", "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, "R.data", sep = ""
      )
      #Attention: data has to be generated!
      #load(here(paste0("data/", nameLoad))) 
      #has a large array "estArrayMean"
      #Update Tibble with empirical variance of sliding based estimator for m
      cat(paste("Sliding est tibble: Länge m: ",length(mVec), "; Länge y: ", 
                length(estArrayMean[,,1,1] %>% apply(c(2), var)), "\n"))
      updTipSlidEst <- tibble(
        gamma = gamma, tsModel = ts_mod, marginal = distr,
        m = mVec, r = rVec[1],
        k = "slidingEst",
        y = estArrayMean[,,1,1] %>% apply(c(2), var)
      )
      cat(paste("Bst tibble: Länge m: ",length(mVec)*length(kVec), "; Länge y: ", 
                length((estArrayMean[,,,seq(2,B+1)] %>% 
                          apply(c(1,2,3), var)) %>% 
                         apply(c(2,3), mean) %>% as.numeric()), "\n"))
      updTibBst <- tibble(
        gamma = gamma, tsModel = ts_mod, marginal = distr, 
        m = rep(mVec, length(kVec)), r = rVec[1],
        k = as.factor(rep(kVec, each = length(mVec))),
        y = (estArrayMean[,,,seq(2,B+1)] %>% 
               apply(c(1,2,3), var)) %>% apply(c(2,3), mean) %>% as.numeric()
      )
      bstTib <- bind_rows(bstTib,updTipSlidEst, updTibBst)
    }
  }
}
nameSave <- "bstMeanEstTib"
save(bstTib, 
     file = here(paste0("data/", nameSave)))