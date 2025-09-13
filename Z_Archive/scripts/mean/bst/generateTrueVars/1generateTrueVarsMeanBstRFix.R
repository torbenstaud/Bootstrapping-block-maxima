library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(tidyverse)

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parametersMeanBstRFixTrueVars.R



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
n <- m * r
cat("Parameter gamma, ts_mod, distr, beta, m, r= ", 
    gamma, ts_mod, distr, beta, m,r, "\n \n", 
    n, m, r, distr, ts_mod, gamma, 10^3, beta)


t0 <- Sys.time()

vars <- trueMeanEstVarC(
  m, r, distr, ts_mod, gamma, 10^6, beta
)

tdel <- difftime(Sys.time(), t0, units = "mins")
cat(paste("It took for:", gamma,ts_mod,distr,beta,m,r, tdel,
          " minutes."))
nameSave <- paste0("trueVar",gamma,ts_mod,distr,beta,m,r)


save(vars, file = here(paste0("data/", nameSave)))




cat("\n r was fixed as n = ",r)