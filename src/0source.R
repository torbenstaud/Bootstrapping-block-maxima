#kmax Funktion: benötigt Rcpp Funktion aus "0sourceC.cpp"----
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


#Eine Table Funktion, die Matrix ausgibt----
#einfach table Funktion in Rcpp
tableC <-
  "#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector tableC(NumericVector v){
  return(table(v));
}"
sourceCpp(code = tableC)
#table als Matrix
#benötigt die oben gesourced Rcpp Funktion "tableC"
cTableMat <- function(x){
  tab <- tableC(x)
  #Liste speichern
  return(
    matrix(data = c(as.double(names(tab)), tab), nrow = 2, byrow = T)
  )
}

#table function erstellen, die eine Liste mit Matrizen (tables) ausgibt----
#: Input.  1: Numeric Vector xx: n Daten, 2: integer l: Größe der Teilblöcke
# Output. Eine Liste mit Matrizen, die die Tables für die jeweiligen disjunkten 
#   Blöcke an l Beobachtungen beinhalten.
# Achtung: Der letzte Block kann eventuell weniger als l Beobachtungen beinhalten,
#   falls l nicht n teilt.

lTable <- function(x, l, n = 0){
  if(n == 0){
    n <- length(x)
  }
  mL <- ceiling(n/l)
  resList <- vector(mode = "list", length = mL)
  if(mL > 1){
    for(indI in seq(1,mL - 1)){
      resList[[indI]] <- cTableMat(x[seq((indI -1)*l+1, indI*l)])
    }
  }
  resList[[mL]] <- cTableMat(x[seq((mL-1)*l+1, n)])
  return(resList)
}





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

#fehlen noch Tests:
#mle (alpha, sigma) aus einem numerischen vektor berechnen; nicht für bootstrapping geeignet
mleFre <- function(xx){
  critFun <- function(a){
    psiKC(a, xx)
  }
  hatA <- uniroot(critFun, c(0.01, 50))$root #lower und upper aus Axels Skript übernommen
  return(c(hatA, hatSigmaCPlain(hatA, xx)))
}


#mle shape aus einem numerischen vektor berechnen; nicht für bootstrapping geeignet

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
    if(indN %% 10000 == 0){
      t1 <- Sys.time()
      timeDiff <- difftime(t1, t0, units = "mins")
      cat("\n \n Repetition: ", indN, " out of ", N,
          "\n Executing time: ", timeDiff, " Minutes.", 
          "\n Remaining: ", (N/indN-1)*timeDiff," Minutes.",
          file = here(nameProgress), append = T)
    }
  }
  return(var(ests))
}



#MLE aus einem ltable berechnen----
#Input: ltable (also Liste mit tables von Unterblöcken)
#benötigt Rcpp Funktionen "psiKCTab" und "hatSigmaC" aus Datei "CFcnts.cpp"

mleFreTab <- function(xx){
  vecX <- unlist(xx)
  critFun <- function(a){
    psiKCTab(a, vecX)
  }
  hatA <- uniroot(critFun, c(0.01, 50))$root #lower und upper aus Axels Skript übernommen
  return(c(hatA, hatSigmaC(hatA, vecX)))
}

#selbe Funktion aber basierend auf der Vektorversion der tables
mleFreTabVec <- function(xx){
  vecX <- unlist(xx)
  critFun <- function(a){
    psiKCVec(a, vecX)
  }
  hatA <- uniroot(critFun, c(0.01, 50))$root #lower und upper aus Axels Skript übernommen
  return(c(hatA, hatSigmaCVec(hatA, vecX)))
}

#Resampling Bootstrap, der ein Replicate zurückgibt----
#Er funktioniert mit allen Listenstrukturen; d.h.: Liste mit Einträgen, die 
# gebootstrapped werden sollen (also zB tableVecList, tableMatrixList)
kBootstrap <- function(xx){
  len <- length(xx)
  inds <- sample(seq(1,len), replace = T)
  bootSamp <- vector(mode = "list", length = len)
  for(i in seq(1, len)){
    bootSamp[[i]] <- xx[[inds[i]]]
  }
  return(bootSamp)
}
