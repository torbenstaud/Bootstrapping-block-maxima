########################################################
# Run file - Generate Bootstrapped mean of block maxima statistics
########################################################
#!/usr/bin/env Rscript
#necessary libraries
library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths 
library(tidyverse) #tibbles

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parameters.R in the same folder as this R.file

r <- rVec[1] #fixed r for this study
#del
gammaVec <- c(-0.2, -0.1, 0, 0.1, 0.2) #frechet Shape in the limit
tsModVec <- c(3) # 1 =iid, 2 = m-dep, 3 = armax
distrVec <- c(1,5) #1= GPD-gamma , 2= Normal, 3=frechet, 4 = pareto, 5 = GEV
betaVec <- c(0, 0.25, 0.5, 0.75) # 0 = iid
#del
fullTrueVarMeanTib <- tibble()
for(gamma in gammaVec){
  for(ts_mod in tsModVec){
    for(distr in distrVec){
      for(beta in betaVec){
        for(m in mVec){
          nameLoad <- paste0("trueTmpVarsTib",gamma,ts_mod,distr, beta, m,r)
          load(here(paste0("data/r", r, "/", nameLoad)))
          updTib <- trueTmpVarsTib # %>% mutate(beta = beta)
          fullTrueVarMeanTib <- bind_rows(
            fullTrueVarMeanTib,
            updTib)
        }
      }
    }
  }
}

save(fullTrueVarMeanTib, file = here(paste0("data/r", r, "/fullTrueVarMeanTib")))


