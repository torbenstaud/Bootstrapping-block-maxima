########################################################
# Run file - Generate Bootstrapped Frechet-MLE statistics
########################################################
#!/usr/bin/env Rscript
#necessary libraries
library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parametersRFix.R in the folder of the file

#Initialize HPC index from for loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#load different parameters and create grid for the HPC to loop over----
combs <- expand.grid(gammaVec, tsModVec, distrVec, betaVec)
# initialize the parameters from the HPC loop ----
gamma <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
cat("Parameter gamma, ts_mod, distr, beta = ", gamma, ts_mod, distr, beta)
r <- rVec[1] #fixed r for this study


#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
nameProgress <- paste(
  "Progress_MleEst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr,
  "_beta",beta, ".txt", sep = ""
)


nameSave <- paste(
  "MleEst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
  "_beta", beta, ".Rdata", sep = ""
)



dims <- c(N, length(mVec), length(kVec))


#estimator specific
#Which estimator?
estFun <- meanCTabVec
#parameter arrays
meanArray <- array(dim = dims)


#Generating loop
t0 <- Sys.time()
if(TRUE){
  for(NSeed in seq(1,N)){
    for(indM in seq_along(mVec)){
      n <- nVec[indM]
      obs <- pmax(rtsC(n, r, distr, ts_mod, gamma, NSeed, beta), .Machine$double.eps)
      for(indK in seq_along(kVec)){
        if(kVec[indK] == 0){#simple sliding Bootstrap (inkonsistent) Fall
          kMaxSamp <- 
            lTableVec(kMaxC(obs, r, kVec[indK]), l = r)
        }else if(kVec[indK] == 1){#disjoint blocks case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, r, kVec[indK]), l = 1)
        }else{#true kMax case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, r, kVec[indK]), l = r*kVec[indK])
        }
        
        est <- estFun(unlist(kMaxSamp))
        
        meanArray[NSeed, indM, indK] <- est
      }
    }
    if(NSeed %% 100 == 0){
      t1 <- Sys.time()
      timeDiff <- difftime(t1, t0, units = "mins")
      cat("\n \n Repetition: ", NSeed, " out of ", N,
          "\n Executing time: ", timeDiff, " Minutes.", 
          "\n Remaining: ", (N/NSeed-1)*timeDiff," Minutes.",
          file = here(nameProgress), append = T)
    }
  }
}
save(meanArray, file = here(paste0("data/", nameSave)))