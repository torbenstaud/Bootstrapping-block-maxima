//CFctns.cpp
#include <Rcpp.h>
using namespace Rcpp;
#include <cstdlib>  // Für strtod
//Abschnitt 0: Allgemine Hilfsfunktionen
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector convert_to_doubles(CharacterVector char_vec) {
  int n = char_vec.size();
  NumericVector num_vec(n);
  
  for(int i = 0; i < n; ++i) {
    std::string str = Rcpp::as<std::string>(char_vec[i]);
    char* endptr;
    double val = std::strtod(str.c_str(), &endptr);
    
    // Überprüfen, ob die Konvertierung erfolgreich war
    if (endptr == str.c_str()) {
      // Keine Konvertierung möglich, ungültige Zahl
      num_vec[i] = NA_REAL;
    } else {
      num_vec[i] = val;
    }
  }
  
  return num_vec;
}

// ABSCHNITT 1: HILFSFUNKTIONEN UM ZEITREIHEN ZU SAMPLEN
// set seed
// sets the r-seed
// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}


//' the quantile function of the GPD-gamma distr. Takes double as input!
 // [[Rcpp::export]]
 double qgpd(double q, double shape){
   double z;
   if(shape == 0){
     z = - std::log(1-q);
   }else{
     z = ( std::pow(1-q,-shape) - 1)/shape;
   }
   return(z);
 }
 
 /* vectorized version of the quantile function of the GPD gamma distr. 
  Takes NumericVector as input and also output
  */
 // [[Rcpp::export]]
 Rcpp::NumericVector qgpdVec(Rcpp::NumericVector q, double shape){
   Rcpp::NumericVector z (q.length());
   if(shape == 0){
     for(int ind = 0; ind < q.length(); ind++){
       z[ind] = - std::log(1-q[ind]);
     }
   }else{
     for(int ind = 0; ind < q.length(); ind++){
       z[ind] = ( std::pow(1-q[ind],-shape) - 1)/shape;
     }
   }
   return(z);
 }
// [[Rcpp::export]]
double qFreSc(double q, double shape){
  double x = std::pow(-std::log(q), -1/shape);
  return(x);
}

// [[Rcpp::export]]
Rcpp::NumericVector qFre(Rcpp::NumericVector q, double shape){
  int len = q.length();
  Rcpp::NumericVector xx (len);
  for(int ind = 0; ind < len; ind++){
    xx[ind] = qFreSc(q[ind], shape);
  }
  return(xx);
}

// [[Rcpp::export]]
Rcpp::NumericVector pFre(Rcpp::NumericVector x, double shape){
  int len = x.length();
  Rcpp::NumericVector qq (len);
  for(int ind = 0; ind < len; ind++){
    qq[ind] = std::exp( -std::pow(1/x[ind],shape));
  }
  return(qq);
}
// [[Rcpp::export]]
double qParSc(double q, double shape){
  double x = std::pow(1 - q, -1/shape);
  return(x);
}
// [[Rcpp::export]]
Rcpp::NumericVector qPar(Rcpp::NumericVector q, double shape){
  int len = q.length();
  Rcpp::NumericVector xx (len);
  for(int ind = 0; ind < len; ind++){
    xx[ind] = qParSc(q[ind], shape);
  }
  return(xx);
}

/*scalar quantile function of GEV-gamma*/
// [[Rcpp::export]]
double qGevSc(double q, double shape){
  double x  = 0;
  if(shape == 0){
    x = -std::log(-std::log(q));
  }else{
    x = (std::exp(-shape*std::log(-std::log(q))) - 1)/shape;
  }
  return(x);
}

/*Vector valued quantile function of GEV-gamma*/
// [[Rcpp::export]]
Rcpp::NumericVector qGev(Rcpp::NumericVector q, double shape){
  int len = q.length();
  Rcpp::NumericVector xx (len);
  for(int ind = 0; ind < len; ind++){
    xx[ind] = qGevSc(q[ind], shape);
  }
  return(xx);
}


// [[Rcpp::export]]
Rcpp::NumericVector rArmaxC(int n, double beta, int marginal, double shape, int burnin, 
                            int seed){
  set_seed(seed);
  Rcpp::NumericVector datFre = qFre(Rcpp::runif(n+burnin), 1);
  Rcpp::NumericVector xx (n+burnin);
  xx[0] = datFre[0];
  for(int ind = 0; ind < n+burnin-1; ind++){
    xx[ind+1] = std::max(beta * xx[ind], (1-beta) * datFre[ind+1]);
    //if(xx[ind+1] < std::pow(10, -7)){
    //  //std::cout << "Wert klein" << xx[ind+1];
    //  xx[ind+1] = 0;
    //}else if(xx[ind+1] > std::pow(10, 3)){
    //  //std::cout << "Wert groß" << xx[ind+1];
    //  xx[ind+1] = max(xx[Rcpp::seq(0, ind)]);
    //}
  }
  Rcpp::NumericVector resVec = xx[Rcpp::seq(burnin, n + burnin-1 )];
  if(marginal == 1){ // GPD Case
    resVec = qgpdVec(pFre(resVec, 1), shape);
  }else if(marginal == 2){// Norm case
    resVec = Rcpp::qnorm(pFre(resVec, 1));
  }else if(marginal == 3){// Frechet case
    resVec = qFre(pFre(resVec, 1), shape);
  }else if(marginal == 4){// Pareto case
    resVec = qPar(pFre(resVec, 1), shape);
  }else if (marginal == 5){// GEV case
    resVec = qGev(pFre(resVec, 1), shape);
  }
  return(resVec);
}
// [[Rcpp::export]]
// Sample a certain m-dependent time series based on a sliding sum window

Rcpp::NumericVector rMDepSumC(
    int n, int m, int marginal, double gamma, int seed
){
  Rcpp::NumericVector xx (n);
  set_seed(seed);
// m-dependent time series
    Rcpp::NumericVector eps = Rcpp::rnorm(n+m-1, 0.0, std::sqrt(1/(double)m) );
    for(int i = 0; i < n; i++){
      xx[i] = Rcpp::sum(eps[Rcpp::seq(i, i + m -1)]); //now N 0,1 distr
      if( marginal == 1){// GPD-gamma case
          xx[i] = qgpd(R::pnorm(xx[i], 0, 1, true, false ), gamma);
      }else if(marginal == 3){ //Frechet case
          xx[i] = qFreSc(R::pnorm(xx[i], 0, 1, true, false ), gamma);
      }else if(marginal == 4){ //Pareto case
        xx[i] = qParSc(R::pnorm(xx[i], 0, 1, true, false ), gamma);
      }
      //else if(marginal == 2){ // NORM-case
        // Nothing to do here, since xx i is already N 0,1 distributed
      //}
    }
  return(xx);
}

 // creates a time series sample from the following models: 1: i.i.d.
 //   2: m-dependent time series, 3: beta-armax time series with GPD-gamma or 
 //   NORM marginals.
 // manual: marginal = 1 means: GPD gamma marginal; marginal = 2 means: NORM marginal
 //  tsMod: 1= IID, 2= MDEP-m, 3=ARMAX-beta [missing!]
// [[Rcpp::export]]
Rcpp::NumericVector rtsC(
    int n, int m, int marginal, int tsMod, double gamma, 
    int seed, double beta
){
  Rcpp::NumericVector xx (n);
  if(tsMod == 1){ // iid observations
    set_seed(seed);
    if(marginal == 1){ // GPD-gamma case
      xx = qgpdVec( Rcpp::runif(n, 0, 1), gamma);
    }else if(marginal == 2){ // NORM case
      xx = Rcpp::rnorm(n, 0, 1);
    }else if(marginal == 3){ // Frechet case
      xx = qFre(Rcpp::runif(n,0,1), gamma);
    }else if(marginal == 4){ // Pareto case
      xx = qPar(Rcpp::runif(n,0,1), gamma);
    }else if(marginal == 5){ // GEV case
      xx = qGev(Rcpp::runif(n,0,1), gamma);
    }
  }
  else if(tsMod == 2){ // m-dependent time series
    xx = rMDepSumC(n, m, marginal, gamma, seed);
  }
  else if (tsMod == 3){ // ARMAX-beta time series
    xx = rArmaxC(n,beta,marginal, gamma, 100, seed);
    }
  return(xx);
}



 // gives the "true" mean of a block maximum of a MDep time series excerpt with 
 //  GPD-gamma marginals.
 // missing: specify tsMod argument: "IID" or "MDEP" (later on maybe "ARMAX")
// [[Rcpp::export]]
double getTruthC(int r, int marginal, int tsMod, 
                 double gamma, int N, double beta){
  Rcpp::NumericVector bms (N);
  for(int ind = 0; ind < N; ind++){
    bms[ind] = Rcpp::max(rtsC(r, r, marginal, tsMod, gamma, ind, beta));
    //bms[ind] = Rcpp::max(rMDepSumC(r, r, marginal, tsMod, gamma, ind));
  }
  return(Rcpp::mean(bms));
}
// [[Rcpp::export]]
// calculate the "true" variance of a block maximum from the specified time series
double getTrueVarC(int r, int marginal, int tsMod, double gamma, int N, double beta){
  double mu = getTruthC(r, marginal, tsMod, gamma, N, beta);
  double v = 0;
  Rcpp::NumericVector bms (N);
  for(int ind = 0; ind < N; ind ++){
    v += std::pow(Rcpp::max(rtsC(r, r, marginal, tsMod, gamma, ind, beta)) - mu, 2);
  }
  return(v/N);
}
// ABSCHNITT 2: HILFSFUNKTIONEN FUER MLE FUNKTION 
//[[Rcpp::export]]
// berechnet psiK, wobei der Vektor xx ein einfacher numerischer Vektor ist
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
}
//[[Rcpp::export]]
// berechnet psiK, wobei der Vektor xx ein unlisted ltable ist
double psiKCTab(double a, Rcpp::NumericVector xx){
  double s1, s2, s3, n;
  s1 = 0;
  s2 = 0;
  s3 = 0;
  n = 0;
  for(int ind = 0; ind < xx.length()/2; ind++){
    s1 += xx[2*ind+1]*std::pow(xx[2*ind], -a)*std::log(xx[2*ind]);
    //std::cout << "s1" << s1;
    s2 += xx[2*ind+1]*std::pow(xx[2*ind], -a);
    s3 += xx[2*ind+1]*std::log(xx[2*ind]);
    n += xx[2*ind+1];
  }
  double res = 1/a + s1/s2 - s3/n;
  return(res);
}
//[[Rcpp::export]]
// berechnet psiK, wobei der Vektor xx ein unlisted ltable ist, 
//  was named vectors beinhaltet
double psiKCVec(double a, Rcpp::NumericVector xx){
  double s1, s2, s3, n;
  s1 = 0;
  s2 = 0;
  s3 = 0;
  n = 0;
  NumericVector vals = convert_to_doubles(xx.names());
  for(int ind = 0; ind < xx.length(); ind++){
    s1 += xx[ind]*std::pow(vals[ind], -a)*std::log(vals[ind]);
    //std::cout << "s1" << s1;
    s2 += xx[ind]*std::pow(vals[ind], -a);
    s3 += xx[ind]*std::log(vals[ind]);
    n += xx[ind];
  }
  double res = 1/a + s1/s2 - s3/n;
  return(res);
}

//[[Rcpp::export]]
// Berechnet hatsigma, wenn hatAlpha bereits berechnet wurde: xx ist ein numerischer vektor -nicht ge lTabled!-
double hatSigmaCPlain(double hatA, Rcpp::NumericVector xx){
  double s, n;
  s = 0;
  n = 0;
  for(int ind = 0; ind < xx.length(); ind++){
    s += std::pow(xx[ind], -hatA);
    n += 1;
  }
  return(std::pow(s/n, -1/hatA));
}


//[[Rcpp::export]]
// Berechnet hatsigma, wenn hatAlpha bereits berechnet wurde: xx ist ein unlisted lTable
double hatSigmaC(double hatA, Rcpp::NumericVector xx){
  double s, n;
  s = 0;
  n = 0;
  for(int ind = 0; ind < xx.length()/2; ind++){
    s += xx[2*ind+1]*std::pow(xx[2*ind], -hatA);
    n += xx[2*ind+1];
  }
  return(std::pow(s/n, -1/hatA));
}

//[[Rcpp::export]]
// Berechnet hatsigma, wenn hatAlpha bereits berechnet wurde: xx ist ein unlisted lTable
double hatSigmaCVec(double hatA, Rcpp::NumericVector xx){
  double s, n;
  s = 0;
  n = 0;
  NumericVector vals = convert_to_doubles(xx.names());
  for(int ind = 0; ind < xx.length(); ind++){
    s += xx[ind]*std::pow(vals[ind], -hatA);
    n += xx[ind];
  }
  return(std::pow(s/n, -1/hatA));
}



// ABSCHNITT 3: kmax Hifsfunktionen
//function_kmaxC.cpp
#include <Rcpp.h>
using namespace Rcpp;
#include "math.h"

// [[Rcpp::export]]
NumericVector cC(NumericVector x, NumericVector y){
  int lenX = x.length(); int lenY = y.length();
  NumericVector res (lenX + lenY);
  for(int ind = 0; ind < lenX; ind++){
    res[ind] = x[ind];
  }
  for(int ind = 0; ind < lenY; ind++){
    res[ind + lenX] = y[ind];
  }
  return(res);
}

// [[Rcpp::export]]
NumericVector seqC(int a, int b){
  // kein Runtime Vorteil
  NumericVector res (b-a+1);
  for(int ind = 0; ind < b-a +1; ind++){
    res[ind] = ind +a;
  }
  return(res);
}


// [[Rcpp::export]]
NumericVector repC(double val, int times){
  NumericVector res (times);
  for(int ind = 0; ind < times; ind++){
    res[ind] = val;
  }
  return(res);
}


// [[Rcpp::export]]
double maxVecC(NumericVector xx){
  //braucht länger :-(
  int len = xx.length();
  double m = xx[0];
  if(len == 1){
    return m;
  }else{
    for(int ind = 1; ind < len; ind++){
      double comp = xx[ind];
      if(m < comp){
        m = comp;
      }
    }
  }
  return m;
}

// [[Rcpp::export]]
NumericVector slidMaxC(
    NumericVector xx, int r
){
  int len = xx.length();
  NumericVector xxL (len + r -1);
  xxL[seqC(0, len -1)] = xx;
  xxL[seqC(len, len + r -2)] = xx[seqC(0, r-2)];
  NumericVector bms (len);
  for(int ind = 0; ind < len; ind++){
    bms[ind] = max(xxL[seq(ind, ind + r - 1)]);
  }
  return(bms);
}

// wie slidMaxC aber ohne Looping der ersten Beobachtungen
// [[Rcpp::export]]
NumericVector slidMaxCNoLoop(
    NumericVector xx, int r
){
  int len = xx.length();
  NumericVector bms (len-r+1);
  for(int ind = 0; ind < len-r+1; ind++){
    bms[ind] = max(xx[seq(ind, ind + r - 1)]);
  }
  return(bms);
}


// [[Rcpp::export]]
NumericVector kMaxTrC(
    NumericVector sample, int r, int k

){
  // k ist hier dann schon echt zwischen 1 und m; also true kMax Fall!
  int n; 
  double m, mk;
  n = sample.length();
  m = std::floor(1.0*n/r); 
  mk = std::floor(1.0*n/(k*r));
  NumericVector bms (n);
  for(int ind = 0; ind < mk; ind++){
    NumericVector iBlockInd = seqC((ind)*k*r, (ind+1)*k*r-1);
    bms[iBlockInd] = slidMaxC(sample[iBlockInd], r);
  } 
  //now treat the last block and distinguish two cases, which need different handling:
  //size of last block
  int lbSize = n-k*r*mk;
  //if n is not divisible k*r: 
  if(lbSize != 0){
    //last block indices
    NumericVector lbInd = seqC(mk*k*r, n-1);
    //case 1: size of last block is smaller then r: Then take bm of last block
    // and repeat it last-block size- times
    if(lbSize < r){
      NumericVector lastBlock = sample[lbInd];
      double maxLblock = max(lastBlock);
      bms[lbInd] = repC(maxLblock, lbSize);
    }else{ //size of last block >= r: use generic sliding max on last block
      bms[lbInd] = slidMaxC(sample[lbInd], r);
    }
  }
  return(bms);
}

// ABSCHNITT 4: Funktionen fuer die empirische Varianz
// calculate empirical variance of a unlisted lTable
// [[Rcpp::export]]
double varCTab( Rcpp::NumericVector xx){
  double mBar, v, n;
  int l = xx.length();
  mBar = 0;
  v = 0;
  n = 0;
  for(int ind = 0; ind < l/2; ind++){
    //the even indices correspond to values of observations
    // the uneven indices correspond to the frequency of the preceeding value
    mBar += xx[2*ind+1] * xx[2*ind];
    n += xx[2*ind+1];
  }
  mBar = mBar/n; //arithmetic mean
  for(int ind = 0; ind < l/2; ind++){
    v += xx[2*ind+1]*std::pow(mBar - xx[2*ind],2);
  }
  v = v/(n-1); //empirical variance
  return(v);
}



// calculate empirical variance of a unlisted lTable with named Vectors inside
// [[Rcpp::export]]
double varCTabVec( Rcpp::NumericVector xx){
  double mBar, v, n;
  int l = xx.length();
  mBar = 0;
  v = 0;
  n = 0;
  NumericVector vals = convert_to_doubles(xx.names());
  for(int ind = 0; ind < l; ind++){
    mBar += xx[ind] * vals[ind];
    n += xx[ind];
  }
  mBar = mBar/n; //arithmetic mean
  for(int ind = 0; ind < l; ind++){
    v += xx[ind]*std::pow(mBar - vals[ind], 2);
  }
  v = v/(n-1); //empirical variance
  return(v);
}

// calculate empirical mean of a unlisted lTable
// [[Rcpp::export]]
double meanCTab( Rcpp::NumericVector xx){
  double mBar, n;
  int l = xx.length();
  mBar = 0;
  n = 0;
  for(int ind = 0; ind < l/2; ind++){
    //the even indices correspond to values of observations
    // the uneven indices correspond to the frequency of the preceeding value
    mBar += xx[2*ind+1] * xx[2*ind];
    n += xx[2*ind+1];
  }
  mBar = mBar/n; //arithmetic mean
  return(mBar);
}
// calculate empirical mean of a unlisted lTable with named Vectors inside
// [[Rcpp::export]]
double meanCTabVec( Rcpp::NumericVector xx){
  double mBar, n;
  int l = xx.length();
  NumericVector vals = convert_to_doubles(xx.names());
  mBar = 0;
  n = 0;
  for(int ind = 0; ind < l; ind++){
    mBar += xx[ind] * vals[ind];
    n += xx[ind];
  }
  mBar = mBar/n; //arithmetic mean
  return(mBar);
}

//calculate true variance of mean bm sliding based estimator
// [[Rcpp::export]]
double trueMeanEstVarC(int m, int r, int distr, int ts_mod, double gamma, int N, double beta){
  NumericVector ests (N);
  for(int indN = 0; indN < N; indN++){
    ests[indN] = Rcpp::mean(slidMaxC(rtsC(m*r, m, distr, ts_mod, gamma, indN, beta), r));
  }
  return(Rcpp::var(ests));
}





