########################################################
# Run file - Process Bootstrapped blockmaxima Frechet MLE
#   for the fixed n (sample size) case
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
# -0parametersMleEstNFix.R in the same folder as this R.file

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
cat("In this experiment n was fixed \n")


#Compare estimation of frechet shape----

#variances for the respective block maxima estimators----
#initialize tibbles for the loop
varTibMleEst <- tibble() #tibble for the estimated alpha/scale

for(alphaInd in seq_along(alphaVec)){
  alpha <- alphaVec[alphaInd]
  trueShape <- alpha
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      for(betaInd in seq_along(betaVec)){
        beta <- betaVec[betaInd]
        nameLoad <- paste(
          "MleEst_Alpha",alpha, "_Ts", ts_mod, "_Distr", distr, 
          "_beta", beta, ".Rdata", sep = ""
        )
        load(here(paste0("data/",nameLoad)))
        vars <- as.numeric(apply(shapeArray, c(2,3), var))
        biases <- as.numeric(apply(shapeArray, c(2,3), mean)) - trueShape
        updTib <-
          tibble(
            alpha = alpha, tsMod = ts_mod, marginal = distr, beta = beta,
            m = rep(mVec, length(kVec)), r = rep(rVec, length(kVec)), 
            k = as.factor(rep(kVec, each = length(mVec))),
            varEst = vars, biasEst = biases
          ) %>% mutate(mse = varEst + biasEst^2)
        varTibMleEst <- bind_rows(varTibMleEst, updTib)
      }
    }
  }
}

save(varTibMleEst,
     file = here(paste0("data/processedDataN",n,".RData")))
cat(paste0("\n We had n = ", n))
