########################################################
# Run file - Generate Bootstrapped Mean statistics (for r fixed)
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
# -0parametersBstRFix.R in the folder of the file

#Initialize HPC index from for loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#load different parameters and create grid for the HPC to loop over----
combs <- expand.grid(gammaVec, tsModVec, distrVec, betaVec, mVec)
# initialize the parameters from the HPC loop ----
gamma <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
m <- combs[index, 5]
cat("Parameter gamma, ts_mod, distr, beta, m = ", gamma, ts_mod, distr, beta, m)
r <- rVec[1] #fixed r for this study

n <- r*m
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
nameProgress <- paste(
  "Progress_MeanBst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr,
  "_beta",beta, "_M", m, ".txt", sep = ""
)


nameSave <- paste(
  "MeanBst_R", r, "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, 
  "_beta", beta, "_M", m, ".Rdata", sep = ""
)



dims <- c(N, length(kVec), B)


#estimator specific
#Which estimator?
estFun <- meanCTabVec
#parameter arrays
estArray <- array(dim = dims)


#Generating loop
t0 <- Sys.time()
if(TRUE){
  for(NSeed in seq(1,N)){
      obs <- rtsC(n, r, distr, ts_mod, gamma, NSeed, beta)
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
        
        
        # Bootstrapping loop
        for(BootSeed in seq(1, B)){
          set.seed(BootSeed)
          kBootSamp <- unlist(kBootstrap(kMaxSamp))
          est <- tryCatch(estFun(kBootSamp), error = function(cond){
            cat("\n \n Loop at gamma = ", gamma, "ts_mod = ", ts_mod, " distr = ",
                distr, " NSeed = ", NSeed,", indM=", indM, "indK = ",indK, 
                "BootSeed = ", BootSeed, "\n Exact Error: \n")
            message(conditionMessage(cond))
            return(c(NA,NA))
          })
          estArray[NSeed, indK, BootSeed] <- est
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
save(estArray, file = here(paste0("data/", nameSave)))