
#!/usr/bin/env Rscript
#necessary libraries
library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(tidyverse)

#sourcing of relevant files----
Rcpp::sourceCpp("/gpfs/project/tosta102/bootstrap/0sourceC.cpp")  #c Functions
source("/gpfs/project/tosta102/bootstrap/0source.R") #R Functions
source(here("0parametersBstRFix.R")) #parameters

bstMeanTib <- tibble()

#loop over hpc param
for(gammaInd in seq_along(gammaVec)){
  gamma <- gammaVec[gammaInd]
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
            "MeanBst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
            "_beta", beta, "_M", m, ".Rdata", sep = ""
          ) #name is estArray
          load(here(paste0("data/",nameLoad)))
          bstMeanTib <- bind_rows(
            bstMeanTib,
            tibble(
              gamma = gamma, tsModel = ts_mod, marginal = distr, beta = beta,
              m = m, r = r, k = kVec, 
              varBst = estArray %>% apply(c(1,2), var) %>% apply(c(2), mean)
            )
          )
          t1 <- Sys.time()
          timeDiff <- difftime(t1, t0, units = "mins")
          cat("\n \n Ein Durchlauf hat gebraucht:", timeDiff, " Minutes. \n
              Insgesamt wird es also ca. ",timeDiff*280, " Minutes brauchen.",
              file = "progressProc.txt", append = T)
        }
        
      }
    }
  }
}
nameSave <- "bstMeanTib"
save(bstMeanTib, file = here(paste0("data/", nameSave)))