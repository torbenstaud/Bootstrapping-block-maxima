########################################################
# Run file - Generate Bootstrapped Frechet-MLE confidence bands 
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
# niveau of CI----
niv <- 0.05


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
ciFreTibNfix <- tibble()
tmpFullArray <- array(dim = c(N, length(kVec), B +1))
#loop over hpc param
for(alphaInd in seq_along(alphaVec)){
  alpha <- alphaVec[alphaInd]
  for(tsInd in seq_along(tsModVec)){
    ts_mod <- tsModVec[tsInd]
    for(distrInd in seq_along(distrVec)){
      distr <- distrVec[distrInd]
      for(betaInd in seq_along(betaVec)){
        beta <- betaVec[betaInd]
        #load the estimation array: name shapeArray
        nameLoadEst <- paste(
          "MleEst_Alpha",alpha, "_Ts", ts_mod, "_Distr", distr, 
          "_beta", beta, ".Rdata", sep = ""
        )
        load(paste0("/gpfs/project/tosta102/bootstrap/mle/scripts/fixedN/est/data/",
                    nameLoadEst))
        shapeArrayEst <- shapeArray
        rm(shapeArray)
        
        for(mInd in seq_along(mVec)){
          
          m <- mVec[mInd]
          r <- rVec[mInd]
          #load the bootstrap array: name shapeArray: dim[N,k,B]
          nameLoadBst <- paste(
            "MleBst_N", n, "Alpha",alpha, "_Ts", ts_mod, "_Distr", distr, 
            "_beta", beta, "_M", m, ".Rdata", sep = ""
          )
          load(here(paste0("data/",nameLoadBst)))
          shapeArrayBst <- shapeArray
          rm(shapeArray)
          tmpFullArray[,,seq(2,B+1)] <- shapeArrayBst
          tmpFullArray[,,1] <- shapeArrayEst[,mInd,]
          updTib <- tibble(
            alpha = alpha, tsModel = ts_mod, marginal = distr,
            beta = beta, m = m, r = r, k = kVec
          ) %>%
            bind_cols(tmpFullArray %>% 
                        apply(c(1,2), function(xx){bootCalcAux(xx, niv)}) %>%
                      ciAna(alpha)
                        )
          ciFreTibNfix <- bind_rows(ciFreTibNfix, updTib)
          
          
        }
        
      }
    }
  }
}

save(ciFreTibNfix,
     file = here(paste0("data/CiFreNfix")))