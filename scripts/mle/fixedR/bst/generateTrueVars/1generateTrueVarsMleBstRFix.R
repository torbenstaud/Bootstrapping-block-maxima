library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(tidyverse)

#sourcing of relevant files----
Rcpp::sourceCpp("/gpfs/project/tosta102/bootstrap/0sourceC.cpp")  #c Functions
source("/gpfs/project/tosta102/bootstrap/0source.R") #R Functions
source(here("0parametersMleBstRFixTrueVars.R")) #parameters

#Initialize HPC index from foor loop----
args <-  commandArgs(trailingOnly=TRUE)
index <- as.numeric(args[1])
#parametrisiere die verschiedenen Kombinationen der 3 Parameter----
combs <- expand.grid(alphaVec, tsModVec, distrVec, betaVec, mVec)
# weise die HPC Parameter dem Index zu----
alpha <- combs[index, 1]; 
ts_mod <- combs[index, 2]; 
distr <- combs[index, 3]
beta <- combs[index, 4]
m <- combs[index, 5]
n <- m * r
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




cat("\n r was fixed as n = ",r)