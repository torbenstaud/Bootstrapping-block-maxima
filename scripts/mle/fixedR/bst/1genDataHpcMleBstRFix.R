########################################################
# Run file - Generate Bootstrapped Frechet-MLE statistics 
#   for the fixed r (block size) case
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
# -0parametersBstRFix.R in the same folder as this R.file

#Initialize HPC index from for loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#load different parameters and create grid for the HPC to loop over----
combs <- expand.grid(alphaVec, tsModVec, distrVec, betaVec, mVec)
# initialize the parameters from the HPC loop ----
alpha <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
m <- combs[index, 5]
cat("Parameter alpha, ts_mod, distr, beta, m = ", alpha, ts_mod, distr, beta, m)
r <- rVec[1] #fixed r for this study

n <- r*m
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
nameProgress <- paste(
  "Progress_MleBst_R", r, "Alpha",alpha, "_Ts", ts_mod, "_Distr", distr,
  "_beta",beta, "_M", m, ".txt", sep = ""
)


nameSave <- paste(
  "MleBst_R", r, "Alpha",alpha, "_Ts", ts_mod, "_Distr", distr, 
  "_beta", beta, "_M", m, ".Rdata", sep = ""
)



dims <- c(N, length(kVec), B)


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
      obs <- pmax(rtsC(n, r, distr, ts_mod, alpha, NSeed, beta), .Machine$double.eps)
      for(indK in seq_along(kVec)){
        #cat("indK = ", indK)
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
        
        
        # Bootstrapping Schleife
        for(BootSeed in seq(1, B)){
          #cat("BootSeed = ", BootSeed)
          set.seed(BootSeed)
          kBootSamp <- kBootstrap(kMaxSamp)
          est <- tryCatch(estFun(kBootSamp), error = function(cond){
            cat("\n \n Schleife bei alpha = ", alpha, "ts_mod = ", ts_mod, " distr = ",
                distr, " NSeed = ", NSeed,", indM=", indM, "indK = ",indK, 
                "BootSeed = ", BootSeed, "\n Genauer Fehler: \n")
            cat("Fehler! \n Genauer Fehler: \n")
            message(conditionMessage(cond))
            return(c(NA,NA))
          })
          shapeArray[NSeed, indK, BootSeed] <- est[1]
          scaleArray[NSeed, indK, BootSeed] <- est[2]
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
save(shapeArray, scaleArray, file = here(paste0("data/r",r,"/", nameSave)))