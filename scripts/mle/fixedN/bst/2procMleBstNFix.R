########################################################
# Run file - Process Bootstrapped Frechet-MLE statistics 
#   for the fixed n (sample size) case
########################################################

#!/usr/bin/env Rscript
#necessary libraries
library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(tidyverse)

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parametersMleBstNFix.R in the same folder as this R.file

#ATTENTION: 
# Files needed are not provided in the repository. You need to generate them 
#   with the generate* scripts.

bstFreVarTib <- tibble()

#loop over hpc param
for(alphaInd in seq_along(alphaVec)){
  alpha <- alphaVec[alphaInd]
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      for(betaInd in seq_along(betaVec[-4])){ # because of missing last one
        beta <- betaVec[betaInd]
        for(mInd in seq_along(mVec)[-1]){ #because of missing last one
          t0 <- Sys.time()
          m <- mVec[mInd]
          r <- rVec[mInd]
          nameLoad <- paste(
            "MleBst_N", n, "Alpha",alpha, "_Ts", ts_mod, "_Distr", distr, 
            "_beta", beta, "_M", m, ".Rdata", sep = ""
          )
          load(here(paste0("data/",nameLoad)))
          bstFreVarTib <- bind_rows(
            bstFreVarTib,
            tibble(
              alpha = alpha, tsModel = ts_mod, marginal = distr, beta = beta,
              m = m, r = r, k = kVec, 
              varBst = shapeArray %>% apply(c(1,2), var) %>% apply(c(2), mean)
            )
          )
          t1 <- Sys.time()
          timeDiff <- difftime(t1, t0, units = "mins")
          cat("\n \n A run took:", timeDiff, " Minutes. \n
              Total run time approx ",timeDiff*168, " Minutes.",
              file = "2progressProc.txt", append = T)
        }
        
      }
    }
  }
}
nameSave <- "bstFreVarTib"
save(bstFreVarTib, file = here(paste0("data/", nameSave)))