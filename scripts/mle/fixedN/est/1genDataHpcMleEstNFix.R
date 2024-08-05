########################################################
# Run file - Generate Bootstrapped Frechet-MLE statistics
#   for the fixed n (sample size) case
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
# -0parametersMleEstNFix.R in the same folder as this R.file


#Initialize HPC index from foor loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#load different parameters and create grid for the HPC to loop over----
combs <- expand.grid(alphaVec, tsModVec, distrVec, betaVec)
# initialize the parameters from the HPC loop ----
alpha <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
cat("Parameter alpha, ts_mod, distr, beta = ", alpha, ts_mod, distr, beta)
cat("\n n was fixed in this experiment: n = ", n)


#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
nameProgress <- paste(
  "Progress_MleEst_Alpha",alpha, "_Ts", ts_mod, "_Distr", distr,
  "_beta",beta, ".txt", sep = ""
)


nameSave <- paste(
  "MleEst_Alpha",alpha, "_Ts", ts_mod, "_Distr", distr, 
  "_beta", beta, ".Rdata", sep = ""
)



dims <- c(N, length(mVec), length(kVec))


#estimator specific
#Which estimator?
estFun <- mleFreTabVec
#parameter arrays
shapeArray <- array(dim = dims)
scaleArray <- array(dim = dims)


#Generating loop
t0 <- Sys.time()
if(TRUE){
  for(NSeed in seq(1,N)){
    #cat("NSeed = ", NSeed,"\n")
    for(indM in seq_along(mVec)){
      r <- rVec[indM]
      obs <- pmax(rtsC(n, r, distr, ts_mod, alpha, NSeed, beta), .Machine$double.eps)
      for(indK in seq_along(kVec)){
        if(kVec[indK] == 0){#simple sliding Bootstrap (inconsistent) case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, r, kVec[indK]), l = r)
        }else if(kVec[indK] == 1){#disjoint blocks case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, r, kVec[indK]), l = 1)
        }else{#true kMax case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, r, kVec[indK]), l = r*kVec[indK])
        }
        
        est <- estFun(kMaxSamp)
        
        shapeArray[NSeed, indM, indK] <- est[1]
        scaleArray[NSeed, indM, indK] <- est[2]
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
save(shapeArray, scaleArray, file = here(paste0("data/", nameSave)))