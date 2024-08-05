########################################################
# Run file - Process  the true variances of the finite sample Fréchet MLE
#   for the fixed n (sample size) case
########################################################
#!/usr/bin/env Rscript

library(tidyverse)
library(here)
#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parametersMleBstNFix.R in the one level above as this R.file
#main loop----


trueVarsTib <- tibble()
for(alpha in alphaVec){
  for(ts_mod in tsModVec){
    for(distr in distrVec){
      for(beta in betaVec){
        for(m in mVec){
          r <- rVec[which(mVec == m)]
          nameLoad <- paste0("trueVar",alpha,ts_mod,distr,beta,m,r)
          cat("\n load name: ", nameLoad)
          load(here(paste0("data/", nameLoad))) #name of numeric variable is "vars"
          trueVarsTib <- bind_rows(
            trueVarsTib,
            tibble(
              alpha = alpha, tsModel = ts_mod, marginal = distr, beta = beta,
              m = m, r = r, 
              trueVar = vars
            )
          )
        }
      }
    }
  }
}
save(trueVarsTib, file = here(paste0("data/trueVarsTib")))