########################################################
# Run file - Generate estimations and bootstrapped replicates 
#   of block maxima statistics for the mean estimation
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
# -0parameters.R in the same folder as this R.file


# gammaVec #frechet Shape in the limit
# tsModVec # 1 =iid, 2 = m-dep, 3 = armax
# distrVec  #1= GPD-gamma , 2= Normal, 3=frechet, 4 = pareto
# 
# N  #number of experiment repitions
# B  #number of bootstrap repetitions
# beta  #armax dependence parameter
# rVec  #value of r (fixed r setting!)
# mVec  #numbers of effective sample sizes
# kVec  #different k-max estimators

#Initialize HPC index from for loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#load different parameters and create grid for the HPC to loop over----
combs <- expand.grid(gammaVec, tsModVec, distrVec)
# initialize the parameters from the HPC loop ----
gamma <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
cat("Parameter gamma, ts_mod, distr= ", gamma, ts_mod, distr)
r <- rVec[1] #fixed r for this study
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
nameProgress <- paste(
  "Progress_mean_", "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, ".txt", sep = ""
)
nameEst <- paste(
  "meanEst_", "gamma",gamma, "_Ts", ts_mod, "_Distr", distr, "R.data", sep = ""
)



#estimator specifics
#Which estimator?
estFun <- function(xx){#argument xx will be a list with tables
  return(meanCTabVec((unlist(xx))))
}
#parameter arrays
dims <- c(N, length(mVec), length(kVec), B+1)
estArrayMean <- array(dim = dims)

#Generating loop
t0 <- Sys.time()
if(TRUE){
  for(NSeed in seq(1,N)){
    for(indM in seq_along(mVec)){
      obs <- rtsC(rVec[1]*mVec[indM], rVec[1], distr, ts_mod, gamma, NSeed, beta)
      for(indK in seq_along(kVec)){
        if(kVec[indK] == 0){#simple sliding Bootstrap (inconsistent) case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, rVec[1], kVec[indK]), l = rVec[1])
        }else if(kVec[indK] == 1){#disjoint blocks case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, rVec[1], kVec[indK]), l = 1)
        }else{#true kMax case
          kMaxSamp <- 
            lTableVec(kMaxC(obs, rVec[1], kVec[indK]), l = rVec[1]*kVec[indK])
        }
        
        est <- estFun(kMaxSamp) #mean estimator
        estArrayMean[NSeed, indM, indK, 1] <- est
        
        # Bootstrapping loop
        for(BootSeed in seq(2, B+1)){
          set.seed(BootSeed)
          kBootSamp <- kBootstrap(kMaxSamp)
          est <- tryCatch(estFun(kBootSamp), error = function(cond){
            cat("\n \n Schleife bei gamma = ", gamma, "ts_mod = ", ts_mod, " distr = ",
                distr, " NSeed = ", NSeed,", indM=", indM, "indK = ",indK, 
                "BootSeed = ", BootSeed, "\n Genauer Fehler: \n")
            cat("Fehler! \n Genauer Fehler: \n")
            message(conditionMessage(cond))
            return(NA)
          })
          estArrayMean[NSeed, indM, indK, BootSeed] <- est
        }
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
cat(paste("\n We had r = ", r, " and gamma, tsMod, distr = ",gamma,
           ts_mod, distr))
save(estArrayMean, file = here(paste0("data/", nameEst)))