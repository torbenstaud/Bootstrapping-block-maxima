# --------------------------------------------------------------------
# Script Name: sourceCase.R
# Purpose: Contains basic functions needed in the case study
# Author: Bücher (Axel), Staud (Torben)
# Date: 2024-08
# --------------------------------------------------------------------



#function from Bücher for asymptotic variance of frechet mle
###################################################
### Asymptotic Variances for Frechet ML-Estimator:
### Input: shape parameter alpha
### Output: list with two entries, 
###   var.sb = Asy.Variance Matrix Sliding Blocks
###   var.db = Asy.Variance Matrix Disjoint Blocks
###################################################
asyvar <- function(alpha) {
  result <- list()
  # Euler-constant:
  gam <- -digamma(1)
  # moments for sliding blocks estimator:
  mo11 <-  2*log(2)*(psigamma(2,0)^2 + pi^2/6 -  psigamma(2,0)*log(2) + log(2)^2/3 ) + 
    psigamma(2,0)* pi^2/6 - 7*zeta(3)/4  - psi(2)^2
  mo22 <- 2*log(2)-1
  mo33 <- 4*log(2) -2
  mo12 <- pi^2/12-log(2)^2+(1-gam)*(2*log(2)-1)
  mo13 <- pi^2/12*(1+psigamma(2))+log(2)^2-2*psigamma(2)*log(2)+psigamma(2)-7/8*zeta(3) 
  mo23 <- pi^2/12+1-2*log(2)
  # Entries for Matrix Sigma_Y:
  ts11 <- 2*mo11; ts22 <- 2*mo22; ts33 <- 2*mo33; ts13 <- 2*mo13
  ts12 <- -2*mo12; ts23 <- -2*mo23
  # Matrix Sigma_Y:
  mat.ts <- matrix(c(ts11, ts12, ts13, ts12, ts22, ts23, ts13, ts23, ts33), nrow=3, ncol=3)
  # Matrix M:
  M.mat <- 6/pi^2*matrix(c(1, 1-gam, -1, gam-1, -(1-gam)^2-pi^2/6, 1-gam), nrow=2, ncol=3, byrow=T)
  
  # Asymptotic variance sliding blocks:
  asyvar.sb <- M.mat %*% mat.ts %*% t(M.mat)
  # Asymptotic variance disjoint blocks (inverse of Fisher information in the Frechet model):
  asyvar.db <- 6/pi^2 * matrix( c(1, gam-1, gam-1, (1-gam)^2+pi^2/6), nrow=2, ncol=2)
  
  result$var.sb <- asyvar.sb %*% diag(c(alpha^2,alpha^(-2)) ) 
  result$var.db <- asyvar.db %*% diag(c(alpha^2,alpha^(-2)) ) 
  return(result)
}

#Asymptotic variance of shape estimation based on the sliding blocks/db mle
v11 <- function(alpha, type = "sl") {
  
  if(type == "sl"){
    return((asyvar(alpha)$var.sb)[1,1])
  }else{
    return((asyvar(alpha)$var.db)[1,1])
  }
}
v11 <- Vectorize(v11)

#Asymptotic variance of scale estimation based on the sliding blocks/db mle

v22 <- function(alpha, type = "sl") {
  if(type == "sl"){
    return((asyvar(alpha)$var.sb)[2,2])
  }else{
    return((asyvar(alpha)$var.db)[2,2])
  }
}
v22 <- Vectorize(v22)


#function from Staud for basic confidence intervals for the Fréchet Pseudo MLE
# based on circmax block maxima.
###################################################
### Circmax conficence intervals:
### Input: xx: sliding block maxima sample, r: block size; k: k block parameter;
###   niv: asymptotic niveau for confidence interval; B: bootstrap replicate
### Output: list with two entries, 
###   quants = Matrix with lower and upper bounds for shape and scale parameters
###   variance = vector with variances for shape and scale bootstrap estimation
###################################################
ciCircmax <- function(xx, r, k, niv = 0.05, B = 10^3){
  xxL <- lTableVec(xx, l = r*k)
  bstSamp <- matrix(nrow = B, ncol = 2) #col1: shape, col2: scale
  for(indB in seq(1,B)){
    set.seed(indB)
    bstSamp[indB,] <- 
      tryCatch(mleFreTabVec(kBootstrap(xxL)), error = function(cond){
        cat("\n \n Erorr at indB = ", indB, "\n exact Error: ")
        message(conditionMessage(cond))
        return(c(NA,NA))}
      )
  }
  bstSamp <- na.omit(bstSamp)
  quants <-
    apply(bstSamp,
          c(2), 
          function(xx) {quantile(xx, c(1-niv/2, niv/2))}
    ) #col1: shape Qs #col2: scale Quantiles
  #auxiliary function to calculate basic boostrap confidence intervals
  auxFun <- function(xx, xMat){
    resMat <- matrix(nrow = 2, ncol = 2)
    for(ind in seq(1,2)){
      resMat[,ind] <- 2*xx[ind] - xMat[,ind]
    }
    return (resMat)
  }
  est <- mleFre(xx)
  return( list(quants = auxFun(est, quants), 
               variance = 
                 apply(bstSamp, c(2), var)))
}
