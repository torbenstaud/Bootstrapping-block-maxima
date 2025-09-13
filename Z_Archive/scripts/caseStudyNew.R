# --------------------------------------------------------------------
# Script Name: caseStudyNew.R
# Purpose: This scripts produces runtime comparison plots for disjoint blocks 
#          vs circmax bootstrapping.
# Author: Torben Staud
# Date: 2024-08
# --------------------------------------------------------------------

# Clear the environment and console
rm(list=ls())  # Clear all objects from the current R environment
cat("\014")    # Clear the console

# ==== SECTION 1: Load required libraries ======================================
# Description: Load all the R packages required for the script.
#_______________________________________________________________________________
library(quantmod) #download financial data from yahoo finance
library(tidyverse) #tidyverse package suite for data manipulation
require(here) #construct paths that are relative to the project's root directory
library(Rcpp) #enable the use of C++ code within R for performance
library(gsl) # package for zeta functions etc.


# ==== SECTION 2: SOURCE EXTERNAL FILES ===========================================
# Description: Source external scripts which might include custom functions
#              or specific project configurations.
#______________________________________________________________________________
source(here("src/0source.r"))
sourceCpp(here("src/0sourceC.cpp"))
source(here("src/sourceCase.R"))
source(here("src/theme.R"))

# ==== SECTION 3: main script ===========================================
# Description: code for the case study data and plots
#______________________________________________________________________________
start_date <- as.Date("1974-01-01")
end_date <- as.Date("2023-12-31")
data0 <- getSymbols("^GSPC", src = "yahoo", from = start_date, 
                    to = end_date, auto.assign = FALSE)
data0_tibble <- as_tibble(data.frame(Date = index(data0), coredata(data0)))
logrets <- data0_tibble %>% mutate(log_return = log(GSPC.Close/lag(GSPC.Close))) %>% 
  drop_na(log_return) %>% mutate(log_wins = pmax(log_return,0), log_loss = -pmin(log_return,0)) %>% 
  select(c(Date, log_wins, log_loss)) 
plogrets <- logrets$log_wins #wins (positive)
nlogrets <- logrets$log_loss #losses (positive)
# Sample Size:
n <- length(logrets$log_wins)
# Dates:
date_all <- as.Date(index(data0))
date_years <- format(date_all, "%Y")
date_quarters <- quarters(date_all)
date_months <- months(date_all)
date_final <- paste(date_years,":",date_quarters, sep="")
date_final_m <- paste(date_years,":",date_months, sep="")  


# Number of quarters:
l.q <- 50*4 #50 years correspond to 50*4 quarters

# Calculate Index where quarter changes (variable qc), and the quarter number (from 1 to l.q=200) 
# each observation belongs to (variable date.quarters2)
qc <- numeric(l.q)
date_quarters2 <- numeric(n)
date_quarters2[1] <- 1; qc[1] <- 1  # Initialize
z <- 1
for(i in seq(2,n)) {
  if(date_quarters[i] != date_quarters[i-1]) {
    z <- z+1
    qc[z] <- i
  }
  date_quarters2[i] <- z
}
# Calculate Quarterly Maxima (disjoint):
plogrets_dbm <- numeric(l.q)
nlogrets_dbm <- numeric(l.q)
for(i in 1:(l.q-1)) {
  plogrets_dbm[i] <- max(plogrets[qc[i]:(qc[i+1]-1)])
  nlogrets_dbm[i] <- max(nlogrets[qc[i]:(qc[i+1]-1)])
}
plogrets_dbm[l.q] <- max(plogrets[qc[l.q]])
nlogrets_dbm[l.q] <- max(nlogrets[qc[l.q]])




#initialize estimation arrays
pEst.db <- matrix(nrow=l.q - 40, ncol=1)
nEst.db <- matrix(nrow=l.q - 40, ncol=1)
pEst.sb <- matrix(nrow=l.q - 40, ncol=1)
nEst.sb <- matrix(nrow=l.q - 40, ncol=1)
pEst.cb <- matrix(nrow=l.q - 40, ncol=1)
nEst.cb <- matrix(nrow=l.q - 40, ncol=1)
pCiUpper <- matrix(nrow=l.q - 40, ncol=1)
pCiLower <- matrix(nrow=l.q - 40, ncol=1)
nCiUpper <- matrix(nrow=l.q - 40, ncol=1)
nCiLower <- matrix(nrow=l.q - 40, ncol=1)
pDbCiUpper <- matrix(nrow=l.q - 40, ncol=1)
pDbCiLower <- matrix(nrow=l.q - 40, ncol=1)
nDbCiUpper <- matrix(nrow=l.q - 40, ncol=1)
nDbCiLower <- matrix(nrow=l.q - 40, ncol=1)



# length of quarter: (62 rougly corresponds to one trading quarter)
r <- 62 
# Estimation loop (rolling window over the quarters)

#extremal probabilities for the negative log returns----

##determine a treshold value: quantile of the block maxima
treshQuants <- quantile(
  kMaxC(nlogrets, r, 1), c(0.45, 0.5, 0.9, 0.95)
) # ~4 and 6 %

treshold <- treshQuants[2] #median
treshInd <- function(xx){
  ifelse(xx <= treshold, 1, 0)
}

#array for bootstrap variances
nbstVars <- array(dim = c(l.q - 40, ncol = 2)) #col1 = dbVar, col2 = cbVar
if(T){ 
  t0 <- Sys.time()
  for(i in 1:(l.q - 40)) {
    ind <- which( date_quarters2 %in% (i:(i+39)) ) #indices for the i-th to i+39-th quarter
    ndata <- nlogrets[ind] #neg ret data from i-th quarter
    # Disjoint Blocks:
    #nEst.db[i,] <- mean(treshInd(nlogrets_dbm[i:(i+39)]))
    nEst.db[i,] <- ndata %>% kMaxC(r=r, k=1) %>% treshInd() %>%  mean()
    nDbQuantsAndVars <- ndata %>% 
      ciCircmaxProb(r = r, k = 2, mthd = "db", tresh = treshold)
    nbstVars[i,1] <- nDbQuantsAndVars[[2]]
    nDbCiBst <- nDbQuantsAndVars[[1]] 
    nDbCiUpper[i] <- nDbCiBst[2] 
    nDbCiLower[i] <- nDbCiBst[1] 
    # Sliding Blocks:
    #pEst.sb[i] <- mean(treshInd(slidMaxCNoLoop(pdata,r)))
    #nEst.sb[i] <- mean(treshInd(slidMaxCNoLoop(ndata,r)))
    nEst.sb[i] <- mean(treshInd(kMaxC(ndata,r, 0)))
    nQuantsAndVars <- ndata %>% 
      ciCircmaxProb(r = r, k = 2, mthd = "cb", tresh = treshold)
    nbstVars[i,2] <- nQuantsAndVars[[2]]
    nCiBst <- nQuantsAndVars[[1]]
    nCiUpper[i] <- nCiBst[2]
    nCiLower[i] <- nCiBst[1] 
    #cb estimator to check
    nEst.cb[i] <- mean(treshInd(kMaxC(ndata,r,2)))
  }
  print(difftime(Sys.time(), t0, units = "min")) #takes approx 0.4 mins
  # save(
  #   nEst.db, nEst.sb, pEst.cb, nEst.cb,
  #   nCiUpper, nCiLower,
  #   nDbCiUpper, nDbCiLower,
  #   
  #   file = here("data/caseStudyDataExtrProb.R")
  # )
}

#quick and dirty comparison
nEst.cb %>% var()
nEst.db %>% var()
nEst.sb %>% var()
nbstVars %>% apply(c(2), mean) %>% 
  (function(xx) {xx[1]/xx[2]})

cbWidths <- nCiUpper - nCiLower
dbWidths <- nDbCiUpper - nDbCiLower
(dbWidths/cbWidths) %>% mean()




if(F){

#mean of the negative log returns blockmaxima ----

#array for bootstrap variances
nbstVars <- array(dim = c(l.q - 40, ncol = 2)) #col1 = dbVar, col2 = cbVar
if(T){ 
  t0 <- Sys.time()
  for(i in 1:(l.q - 40)) {
    ind <- which( date_quarters2 %in% (i:(i+39)) ) #indices for the i-th to i+39-th quarter
    ndata <- nlogrets[ind] #neg ret data from i-th quarter
    # Disjoint Blocks:
    #nEst.db[i,] <- mean(nlogrets_dbm[i:(i+39)])
    nEst.db[i,] <- ndata %>% kMaxC(r=r, k=1) %>% mean()
    nDbQuantsAndVars <- ndata %>% 
      ciCircmaxMean(r = r, k = 2, mthd = "db")
    nbstVars[i,1] <- nDbQuantsAndVars[[2]]
    nDbCiBst <- nDbQuantsAndVars[[1]] 
    nDbCiUpper[i] <- nDbCiBst[2] 
    nDbCiLower[i] <- nDbCiBst[1] 
    # Sliding Blocks:
    #nEst.sb[i] <- mean(treshInd(slidMaxCNoLoop(ndata,r)))
    nEst.sb[i] <- mean(kMaxC(ndata,r, 0))
    nQuantsAndVars <- ndata %>% 
      ciCircmaxMean(r = r, k = 2, mthd = "cb")
    nbstVars[i,2] <- nQuantsAndVars[[2]]
    nCiBst <- nQuantsAndVars[[1]]
    nCiUpper[i] <- nCiBst[2]
    nCiLower[i] <- nCiBst[1] 
    #cb estimator to check
    nEst.cb[i] <- mean(kMaxC(ndata,r,2))
  }
  print(difftime(Sys.time(), t0, units = "min")) #takes approx 0.4 mins
  # save(
  #   nEst.db,nEst.sb,  nEst.cb,
  #   nCiUpper, nCiLower,
  #   nDbCiUpper, nDbCiLower,
  #   
  #   file = here("data/caseStudyDataExtrProb.R")
  # )
}

#quick and dirty comparison
nEst.cb %>% var()
nEst.db %>% var()
nEst.sb %>% var()
nbstVars %>% apply(c(2), mean) %>% 
  (function(xx) {xx[1]/xx[2]})

cbWidths <- nCiUpper - nCiLower
dbWidths <- nDbCiUpper - nDbCiLower
(dbWidths/cbWidths) %>% mean()


}