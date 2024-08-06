require(here)
require(Rcpp)
require(tidyverse)


source(here("src/0source.R"))
sourceCpp(here("src/0sourceC.cpp"))




#' Bootstrap sliding block maxima characteristics
#'
#' @param xx univariate numeric sample 
#' @param k an integer representing the block parameter in the circmax procedure.
#' @param r an integer representing the block block size to take the max of.
#' @param stat a string representing the statistic of interest. The following 
#' values are possible: mean (of a block maximum), variance (of 
#' a block maximum), quantile (of a block maximum)
#' @param q a double representing the quantile
#' @param B an integer representing the number of bootstrap replicates
#' @param niveau a double representing the confidence interval niveau
#' @param normal a boolean representing whether the confidence interval should 
#' be based on normal approximation.
#'
#' @return a list with 3 values: 'replicates' is a numeric vector containing the 
#' B bootstrap replicates. 'ci' is a numeric vector with two entries: the lower 
#' and upper confidence interval bounds. 'variance' a numeric vector containing 
#' the variance of the bootstrap replicates (bootstrapeped sliding blocks 
#' estimator variance).
#' @export
#'
#' @examples
blockMaxBootstrap <- function(
  xx, k = 2, r, stat, q = 0.95, 
  B = 10^3, niveau = 0.05, normal = F  
){
  if(!is.numeric(xx)){
    stop("xx must be numeric vector!")
  }else if(!(floor(k) == k && k > 1 && k < length(xx)/ r)){
    stop("k must be an integer greater than 1 and has to fit r times in xx!")
  }else if(!(floor(r) == r && r > 1 && r <= length(xx)/ k)){
    stop("r has to be an integer greater than 1 and has to fit k times in xx!")
  }else if(!(stat %in% c("mean", "variance", "quantile"))){
    stop("stat has to be 'mean', 'variance' or 'quantile'")
  }else if(!(is.numeric(q) && 0 < q  && q < 1)){
    stop("q has to be a number between 0 and 1!")
  }else if(!(B >= 1 && floor(B) == B)){
    stop("B has to be an integer greater than 0")
  }else if(!(is.numeric(niveau) && 0 < niveau  && niveau < 1)){
    stop("niveau has to be a number between 0 and 1!")
  }else if(!is.logical(normal)){
    stop("normal has to be a boolean!")
  }

  
  res <- list(replicates = NA, ci = NA, variance = NA)
  bstRep <- numeric(B)
  bms <- slidMaxC(xx, r)
  if(stat == "quantile"){
    for(indB in seq_along(1:B)){
      set.seed(indB)
      bstRep[indB] <- 
        quantile(blockBootstrap(bms, k*r), q)
    }
    est <- quantile(bms, q)
  }else{
    
  
    bmsLTab <- lTableVec(bms, l = k*r)
    func <- NA
    if(stat == "mean"){
      func <- weightMeanC
    }else if(stat == "variance"){
      func <- weightVarC
    }
    bmsUnl <- unlist(bmsLTab)
    est <- func(bmsUnl, as.numeric(names(bmsUnl)))
    for(indB in seq_along(1:B)){
      set.seed(indB)
      bstTmp <- unlist(kBootstrapC(bmsLTab))
      bstRep[indB] <- 
        func( bstTmp, as.numeric(names(bstTmp)))
    }
    
  }

  res[["replicates"]] <- bstRep
  res[["variance"]] <- var(bstRep)
  if(normal == F){
    res[["ci"]] <- 2*est - unname(quantile(bstRep, c(1- niveau/2, niveau/2)))
  }else{
    res[["ci"]] <- est - qnorm(c(1-niveau/2, niveau/2), 
                               sd = sqrt(res[["variance"]]))
  }
  return(res)
}



