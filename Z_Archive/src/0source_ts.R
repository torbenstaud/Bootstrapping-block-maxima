########################################################
# source functions for wgen simulation study
########################################################
#!/usr/bin/env Rscript
#necessary libraries
library(Rcpp) #compiling c code
library(evd) #certain distributions
library(here) #relative paths
library(wxgenR) #weather generating
library(magrittr) #pipe operator

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
      resList[[indI]] <- table(x[seq((indI -1)*l+1, indI*l)]) #eig tableC
    }
  }
  resList[[mL]] <- table(x[seq((mL-1)*l+1, n)]) #eig tableC, aber um zu debuggen
  return(resList)
}

##PWM Schätzer (evtl als Anfangswert)
get_pwm <- function(data){
  pwm = compute_pwm(data) #von Rcpp Funktion compute_pwm
  
  g1 <- function(gam){
    ifelse(gam == 0, log(3)/log(2), (3**gam - 1)/(2**gam -1))
  }
  g2 <- function(gam){
    ifelse(gam == 0, 1/log(2), gam/(gamma(1-gam)*(2**gam-1)))
  }
  g3 <- function(gam){
    ifelse(gam == 0, 0.5772156649015328606065120,(1-gamma(1-gam))/gam)
  }
  
  gam = uniroot(function(x) g1(x) - (3*pwm[3] - pwm[1])/(2*pwm[2] - pwm[1]),
                lower = -2, upper = 2)$root
  sigma = g2(gam)*(2*pwm[2] - pwm[1])
  mu = pwm[1] + sigma*g3(gam)
  return(c(mu, sigma, gam))
}
## Abschnitt mit den ML-GEV Funktionen


gev_mle_cpp <- function(data, start = NULL){
  if (!is.numeric(data) || any(!is.finite(data))) {
    stop("Die Eingabedaten müssen numerisch und endlich sein.")
  }
  
  # Startwerte setzen, falls keine angegeben sind
  if(identical(start, c(-100, -100, -100))){
    return(list(
      par = c(-100, -100, -100),
      value = -100,
      method = "error"
    ))
  }
  if (is.null(start)) {
    start_scale <- sqrt(6 * var(data))/pi
    start <- c(mean(data) - 0.58 * start_scale,start_scale,  0)
  }
  
  # Optimierung
  result <- try({optim(
    par = start,
    fn = function(theta) neg_log_likelihood_gev(theta, data), 
    #xx = data,
    method = "BFGS",
    #lower = c(-10, 0.01, -0.9),
    #upper = c(50, 6, 0.9),
    control = list(maxit = 10**5, fnscale = length(data))
  )}, silent = T)
  return(result)
  # # **Fix: Prüfen, ob `result` eine Liste ist und `par` enthält**
  # if (!inherits(result, "try-error") && is.list(result) && "par" %in% names(result)) {
  #   return(list(
  #     par = result$par,
  #     value = result$value
  #   ))
  # }
  # 
  # return(
  # list(
  #   par = c(-100, -100, -100),
  #   value = -100,
  #   method = "error"
  # ))
  
}



#nun den MLE mit der lVec fun bauen:
gev_mle_cpp_lvec <- function(data, start = NULL, method = "L-BFGS-B") {
  # Prüfe Eingabedaten
  if (!is.numeric(data) || any(!is.finite(data))) {
    stop("Die Eingabedaten müssen numerisch und endlich sein.")
  }
  
  # Prüfe Startwerte
  if(identical(start, c(-100, -100, -100))){
    return(list(
      par = c(-100, -100, -100),
      value = -100,
      method = "error"
    ))
  }
  if (is.null(start)) {
    start_scale <- sqrt(6 * varCTabVec(data))/pi
    start <- c(meanCTabVec(data) - 0.58 * start_scale,start_scale,  0)
  } else if (!is.numeric(start) || length(start) != 3 || any(!is.finite(start))) {
    stop("Startwerte müssen ein numerischer Vektor der Länge 3 sein und endlich.")
  }
  
  # Versuch mit "L-BFGS-B"
  resultBfgs <- try({
    optim(
      par = start,
      fn = function(theta) neg_log_likelihood_gev_lvec(theta, data),
      method = "BFGS",
      #lower = c(-10, 0.01, -0.9),
      #upper = c(50, 6, 0.9),
      control = list(maxit = 10**5, fnscale = 10**3) # ohne fnscale probieren, fnscale = sum(data))
    )
  }, silent = TRUE)
  
  if (!inherits(resultBfgs, "try-error")) {
    return(resultBfgs)
  }
  
  
  # Versuch mit "DEoptim"
  cat("Fehler bei 'L-BFGS-B'. Versuche 'DEoptim'.\n")
  resultDe <- try({
    optimRes <- DEoptim(
      fn = neg_log_likelihood_gev_lvec,
      lower = c(-10, 0.01, -0.9),
      upper = c(50, 6, 0.9),
      xx = data
    )$optim
  }, silent = TRUE)
  if (!inherits(resultDe, "try-error")) {
    return(list(
      par = resultDe$bestmem,
      value = resultDe$bestval,
      method = "deOptim"
    ))
  }
  return(
    list(
      par = c(-100, -100, -100),
      value = -100,
      method = "error"
    )
  )
  
  # # Versuch mit "Nelder-Mead"
  # cat("Fehler bei 'L-BFGS-B' und 'deOptim'. Versuche 'Nelder-Mead'.\n")
  # 
  # resultNm <- try({
  #   optim(
  #     par = start,
  #     fn = function(theta) neg_log_likelihood_gev_lvec(theta, data),
  #     method = "Nelder-Mead",
  #     control = list(trace = 0)
  #   )
  # }, silent = TRUE)
  # if (!inherits(resultNm, "try-error")) {
  #   return(list(
  #     par = resultNm$par,
  #     value = resultNm$value,
  #     method = "Nelder-Mead"
  #   ))
  # }
  # 
  # # Fehler bei allen Methoden
  # cat("Fehler bei allen Optimierungsmethoden.")
}



#calculate the T-return level if theta is known (we will use it as plug-in for the GEV MLE)
##arguments: Time(double), theta(GEV parameter vector)
gev_rl <- function(Time, theta){
  if(theta[2] < 0){
    return(-100)
  }
  return(evd::qgev(1- 1/Time, loc = theta[1], scale = theta[2], shape = theta[3]))
}

#Functions for Fréchet MLE----

#berechnet psiK, wobei der Vektor xx ein einfacher numerischer Vektor ist
Rcpp::cppFunction('
double psiKC(double a, Rcpp::NumericVector xx){
  double s1, s2, s3, n;
  s1 = 0;
  s2 = 0;
  s3 = 0;
  n = 0;
  for(int ind = 0; ind < xx.length(); ind++){
    s1 += std::pow(xx[ind], -a)*std::log(xx[ind]);
    //std::cout << "s1" << s1;
    s2 += std::pow(xx[ind], -a);
    s3 += std::log(xx[ind]);
    n += 1;
  }
  double res = 1/a + s1/s2 - s3/n;
  return(res);
}'
)

# Berechnet hatsigma, wenn hatAlpha bereits berechnet wurde: xx ist ein numerischer vektor -nicht ge lTabled!-
Rcpp::cppFunction('double psiKC(double a, Rcpp::NumericVector xx){
  double s1, s2, s3, n;
  s1 = 0;
  s2 = 0;
  s3 = 0;
  n = 0;
  for(int ind = 0; ind < xx.length(); ind++){
    s1 += std::pow(xx[ind], -a)*std::log(xx[ind]);
    //std::cout << "s1" << s1;
    s2 += std::pow(xx[ind], -a);
    s3 += std::log(xx[ind]);
    n += 1;
  }
  double res = 1/a + s1/s2 - s3/n;
  return(res);
}')

#mle (alpha, sigma) aus einem numerischen vektor berechnen; nicht für bootstrapping geeignet
fre_mle_cpp <- function(xx){
  critFun <- function(a){
    psiKC(a, xx)
  }
  hatA <- uniroot(critFun, c(0.01, 50))$root #lower und upper aus Axels Skript übernommen
  return(c(hatSigmaCPlain(hatA, xx), hatA))
}

#calculate the T-return level if theta=(scale, shape) is known (we will use it as plug-in for the Fre MLE)
##arguments: Time(double), theta(scale, shape Fréchet parameter)
fre_rl <- function(Time, theta){
  return(evd::qfrechet(1- 1/Time, loc = 0, scale = theta[1], shape = theta[2]))
}

