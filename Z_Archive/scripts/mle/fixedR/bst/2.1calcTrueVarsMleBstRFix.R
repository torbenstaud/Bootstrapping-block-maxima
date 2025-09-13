
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
# -0parametersBstRFix.R in the same folder as this R.file

#Initialize HPC index from foor loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#parametrisiere die verschiedenen Kombinationen der 3 Parameter----
combs <- expand.grid(alphaVec, tsModVec, distrVec, betaVec, mVec)
# weise die HPC Parameter dem Index zu----
alpha <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
m <- combs[index, 5]
cat("Parameter alpha, ts_mod, distr, beta, m = ", alpha, ts_mod, distr, beta, m)
r <- rVec[1] #fixed r for this study

trueVarsTib <- tibble(
  alpha = alpha, tsModel = ts_mod, marginal = distr,
  m = m, r = r, k = 0,
  varMeanSlidEst = ???
)
#loop over hpc param
for(alphaInd in seq_along(alphaVec)){
  alpha <- alphaVec[alphaInd]
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      for(betaInd in seq_along(betaVec)){
        beta <- betaVec[betaInd]
        for(mInd in seq_along(mVec)){
          t0 <- Sys.time()
          m <- mVec[mInd]
          nameLoad <- paste(
            "MleBst_R", r, "Alpha",alpha, "_Ts", ts_mod, "_Distr", distr, 
            "_beta", beta, "_M", m, ".Rdata", sep = "")
          load(here(paste0("data/r",r,"/",nameLoad)))
          bstFreVarTib <- bind_rows(
            bstFreVarTib,
            tibble(
              alpha = alpha, tsModel = ts_mod, marginal = distr,
              m = m, r = r, k = kVec, 
              varBst = shapeArray %>% apply(c(1,2), var) %>% apply(c(2), mean)
            )
          )
          t1 <- Sys.time()
          timeDiff <- difftime(t1, t0, units = "mins")
          cat("\n \n Ein Durchlauf hat gebraucht:", timeDiff, " Minutes. \n
              Insgesamt wird es also ca. ",timeDiff*168, " Minutes brauchen.",
              file = "progressProc.txt", append = T)
        }
        
      }
    }
  }
}
nameSave <- "bstFreVarTib"
save(bstFreVarTib, file = here(paste0("data/", nameSave)))