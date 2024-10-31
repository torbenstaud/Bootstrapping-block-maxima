
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
# niveau of CI----
niv <- 0.05

############### BIS HIER NEU GEMACHT FÜR MEAN CIs
#auxiliary functions----
##calculate simple bootstrap intervall from B+1 vector with 1st coordinate estimator----
bootCalcAux <- function(xx, niv){
  return(2*xx[1] - quantile(xx[-1], c(1-niv/2, niv/2)))
}

##calculate average ci width from N cis and empirical coverage ----
ciAna <- function(arr, truth){#arr [2,N,kVec] array
  empCov <- 
    apply(arr, c(2,3), function(vec){vec[1] <= truth && vec[2] >= truth}) %>% 
    apply(c(2), mean)
  avgWidth <- apply(arr, c(2,3), function(vec){vec[2] - vec[1]}) %>% 
    apply(c(2), mean)
  return(tibble(empCov = empCov, avgWidth = avgWidth)) #tibble with length kVec
}

#initialize tibble and main loop----
ciMeanTibRfix <- tibble()
tmpFullArray <- array(dim = c(N, length(kVec), B +1))
#loop over hpc param
for(gammaInd in seq_along(gammaVec)){
  gamma <- gammaVec[gammaInd]
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      for(betaInd in seq_along(betaVec)){
        beta <- betaVec[betaInd]
        #load the estimation array: name shapeArray
        nameLoadEst <- paste(
          "MleEst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
          "_beta", beta, ".Rdata", sep = ""
        )
        load(paste0("/gpfs/project/tosta102/bootstrap/mean/scripts/est/data/",
                    nameLoadEst))
        meanArrayEst <- meanArray
        rm(meanArray)
        trueMean <- getTruthC(r, distr, ts_mod, gamma, 10^6, beta)
        t0 <- Sys.time()
        for(mInd in seq_along(mVec)){
          m <- mVec[mInd]
          #load the bootstrap array: name estArray: dim[N,k,B]
          nameLoadBst <- paste(
            "MeanBst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
            "_beta", beta, "_M", m, ".Rdata", sep = ""
          )
          load(here(paste0("data/",nameLoadBst)))
          meanArrayBst <- estArray
          rm(estArray)
          tmpFullArray[,,seq(2,B+1)] <- meanArrayBst
          tmpFullArray[,,1] <- meanArrayEst[,mInd,]
          updTib <- tibble(
            gamma = gamma, tsModel = ts_mod, marginal = distr,
            beta = beta, m = m, r = r, k = kVec
          ) %>%
            bind_cols(tmpFullArray %>% 
                        apply(c(1,2), function(xx){bootCalcAux(xx, niv)}) %>%
                        ciAna(trueMean)
            )
          ciMeanTibRfix <- bind_rows(ciMeanTibRfix, updTib)
          
          
        }
        t1 <- Sys.time()
        timeDiff <- difftime(t1, t0, units = "mins")
        cat("\n \n Ein mVec Durchlauf hat gebraucht:", timeDiff, " Minutes. \n
              Insgesamt wird es also ca. "
            ,timeDiff*length(gammaVec)*length(tsModVec)*length(distrVec)*length(betaVec), 
            " Minutes brauchen.",
            file = "progress2.0Ci.txt", append = T)
        
      }
    }
  }
}

save(ciMeanTibRfix,
     file = here(paste0("data/ciMeanTibRfix")))