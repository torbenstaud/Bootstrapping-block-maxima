########################################################
# Run file - Generate the true variances of the finite sample Fréchet MLE
#   for the fixed n (sample size) case
########################################################
#!/usr/bin/env Rscript

library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(tidyverse)

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parametersMleBstNFixTrueVars.R in the same folder as this R.file

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
r <- rVec[which(mVec == m)]
cat("Parameter alpha, ts_mod, distr, beta, m, r= ", 
    alpha, ts_mod, distr, beta, m,r, "\n \n", 
    n, m, r, distr, ts_mod, alpha, 10^3, beta)


t0 <- Sys.time()

vars <- trueShapeEstVar(
  n = n, m = m, r = r, marginal = distr, tsMod = ts_mod, 
  alpha = alpha, N = 10^6, beta = beta
)

tdel <- difftime(Sys.time(), t0, units = "mins")
cat(paste("It took for:", alpha,ts_mod,distr,beta,m,r, tdel,
          " minutes."))
nameSave <- paste0("trueVar",alpha,ts_mod,distr,beta,m,r)


save(vars, file = here(paste0("data/", nameSave)))




cat("\n n was fixed as n = ",n)