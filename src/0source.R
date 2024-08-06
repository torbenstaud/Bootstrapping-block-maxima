# --------------------------------------------------------------------
# Script Name: sourceCase.R
# Purpose: Contains basic functions needed throughout the study
# Author:Torben Staud
# Date: 2024-08
# --------------------------------------------------------------------



#' Calculate circmax(k) sample of a numeric vector
#'
#' @param sample numeric vector 
#' @param r integer representing the block size
#' @param k integer representing the blocks parameter. k = 0 represents the 
#' sliding blocks case, k = 1 the disjoint case and other values the 
#' circmax (k) cases.
#'
#' @return numeric vector containing the block maxima
#' @export
#'
#' @examples
kMaxC <- function(sample, r, k){
  
  n <- length(sample)
  m <- floor(n/r)
  #Checking preconditions for input is omitted in the sim study
  if(k == 1){ # disjoint bm case
    # if last block not full: fill with NAs
    sample <- c(sample,
                rep(NA, 
                    ifelse(
                      (m +1)*r -n == r , 0,
                      (m +1)*r -n))
    )
    #calculate disjoint blockmaxima and ignore the NAs
    #Note: rcpproll method is (approx 3 times-) faster than apply solution
    bms <- RcppRoll::roll_max(sample, n = r, by = r, na.rm = TRUE)
  }else if (k == 0){#sliding blocks case
    #loop the first r-1 observations to the end to ensure correct weighing
    sample <- append(sample, sample[seq(1, r-1)])
    bms <- RcppRoll::roll_max(sample, n = r)
  }else{# this is the true kMax case
    #k <- as.integer(k)
    bms <- kMaxTrC(sample, r, k)
  }
  return(bms)
}


#table function in Rcpp
tableC <-
  "#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector tableC(NumericVector v){
  return(table(v));
}"
sourceCpp(code = tableC)



#table function erstellen, die eine Liste mit Matrizen (tables) ausgibt----
#: Input.  1: Numeric Vector xx: n Daten, 2: integer l: Größe der Teilblöcke
# Output. Eine Liste mit Matrizen, die die Tables für die jeweiligen disjunkten 
#   Blöcke an l Beobachtungen beinhalten.
# Achtung: Der letzte Block kann eventuell weniger als l Beobachtungen beinhalten,
#   falls l nicht n teilt.


#generate list of tables (named vectors) where the tables are built from l blocks
# of the numeric vector x (often a kMaxed vector)
lTableVec <- function(x, l, n = 0){
  if(n == 0){
    n <- length(x)
  }
  mL <- ceiling(n/l)
  resList <- vector(mode = "list", length = mL)
  if(mL > 1){
    for(indI in seq(1,mL - 1)){
      resList[[indI]] <- tableC(x[seq((indI -1)*l+1, indI*l)])
    }
  }
  resList[[mL]] <- tableC(x[seq((mL-1)*l+1, n)])
  return(resList)
}

#calculate the Fréchet MLE of a numeric vector. returns a vector with two 
# elements: first is the estimated shape, second the estimated scale.
mleFre <- function(xx){
  critFun <- function(a){
    psiKC(a, xx)
  }
  hatA <- uniroot(critFun, c(0.01, 50))$root #lower und upper aus Axels Skript übernommen
  return(c(hatA, hatSigmaCPlain(hatA, xx)))
}


#same as above but without calculating the scale
mleFreShape <- function(xx){
  critFun <- function(a){
    psiKC(a, xx)
  }
  hatA <- uniroot(critFun, c(0.01, 50))$root #lower und upper aus Axels Skript übernommen
  return(hatA)
}

#calculate the true shape estimator variance and return it as a single number
trueShapeEstVar <- function(
    n, m, r, marginal, tsMod, alpha, N, beta){
  ests <- numeric(N)
  #progress 
  nameProgress <- paste(
    "Progress_TrueVar_N", n, "Alpha",alpha, "_Ts", ts_mod, "_Distr", marginal,
    "_beta",beta, "_M", m, "_R", r,  ".txt", sep = ""
  )
  t0 <- Sys.time()
  for(indN in seq(1, N)){
    ests[indN] <- 
      mleFreShape(slidMaxC(rtsC(n, m, marginal, tsMod, alpha, indN, beta), r))
  }
  return(var(ests))
}




#calculate Fréchet MLE parameters via lTables (faster for bootstrapping)
#input: a lTable (list of tables)
#output: parameters as a numeric vector
mleFreTabVec <- function(xx){
  vecX <- unlist(xx)
  critFun <- function(a){
    psiKCVec(a, vecX)
  }
  hatA <- uniroot(critFun, c(0.01, 50))$root #lower und upper aus Axels Skript übernommen
  return(c(hatA, hatSigmaCVec(hatA, vecX)))
}

#Resampling Bootstrap
#Input: lTable (list of tables)
#Output: lTable of same length
kBootstrap <- function(xx){
  len <- length(xx)
  inds <- sample(seq(1,len), replace = T)
  bootSamp <- vector(mode = "list", length = len)
  for(i in seq(1, len)){
    bootSamp[[i]] <- xx[[inds[i]]]
  }
  return(bootSamp)
}


#Resampling Bootstrap based on a numeric vector
#Input: numeric vector
#Output: numeric vector
blockBootstrap <- function(xx, l) {
  n <- length(xx)
  mk <- ceiling(n/l)
  lB <- n - (mk -1)*l
  res <- numeric(mk*l)
  inds <- sample(seq(1,mk), replace = T)
  counter <- 1
  for(ind in seq_along(1:(mk))){
    if(inds[ind] == mk){
      res[seq(counter, counter + lB - 1)] <- 
        xx[seq((mk-1)*l+1, n)]
      counter <- counter + lB
    }else{
      res[seq(counter, counter + l -1)] <-
        xx[seq( (inds[ind] -1)*l+1, inds[ind]*l ) ]
      counter <- counter + l
    }
  }
  return(res[seq(1,n)])
}
