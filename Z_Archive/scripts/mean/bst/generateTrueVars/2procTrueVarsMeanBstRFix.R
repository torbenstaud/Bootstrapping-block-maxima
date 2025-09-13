library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(tidyverse)

#sourcing of relevant files----
# the following files need to be sourced: 
# -0sourceC.cpp root
# -0source.R root
# -0parametersMeanBstRFixTrueVars.R




trueVarsTib <- tibble()
for(gamma in gammaVec){
  for(ts_mod in tsModVec){
    for(distr in distrVec){
      for(beta in betaVec){
        for(m in mVec){
          nameLoad <- paste0("trueVar",gamma,ts_mod,distr,beta,m,r)
          cat("\n load name: ", nameLoad)
          load(here(paste0("data/", nameLoad))) #name of numeric variable is "vars"
          trueVarsTib <- bind_rows(
            trueVarsTib,
            tibble(
              gamma = gamma, tsModel = ts_mod, marginal = distr, beta = beta,
              m = m, r = r, 
              trueVar = vars
            )
          )
        }
      }
    }
  }
}
save(trueVarsTib, file = here(paste0("data/trueVarsTib")))