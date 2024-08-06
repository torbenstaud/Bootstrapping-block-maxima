library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(tidyverse)

#sourcing of relevant files----
Rcpp::sourceCpp("/gpfs/project/tosta102/bootstrap/0sourceC.cpp")  #c Functions
source("/gpfs/project/tosta102/bootstrap/0source.R") #R Functions
source(here("0parametersMleBstRFixTrueVars.R")) #parameters




trueVarsTib <- tibble()
for(alpha in alphaVec){
  for(ts_mod in tsModVec){
    for(distr in distrVec){
      for(beta in betaVec){
        for(m in mVec){
          nameLoad <- paste0("trueVar",alpha,ts_mod,distr,beta,m,r)
          cat("\n load name: ", nameLoad)
          load(here(paste0("data/", nameLoad))) #name of numeric variable is "vars"
          trueVarsTib <- bind_rows(
            trueVarsTib,
            tibble(
              alpha = alpha, tsModel = ts_mod, marginal = distr, beta = beta,
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