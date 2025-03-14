function(xx) {
  n <- length(xx)
  inds <- sample(seq(1,n), replace = T)
  return(xx[inds])
}

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
    method = "L-BFGS-B",
    lower = c(-10, 0.01, -0.9),
    upper = c(50, 6, 0.9),
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
gev_mle_cpp_lvec <- function(data, start = NULL, method = "L-BFGS-B") {
  # Prüfe Eingabedaten
  if (!is.numeric(data) || any(!is.finite(data))) {
    stop("Die Eingabedaten müssen numerisch und endlich sein.")
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
      method = "L-BFGS-B",
      lower = c(-10, 0.01, -0.9),
      upper = c(50, 6, 0.9),
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

ciCircmaxRl <- 
  function(xx, Time = 100, r, k, niv = 0.05, B = 10^3,
           mthd = "cb"){
    #t0 <- Sys.time()
    
    bstSamp <- vector(mode = "numeric", length = B)
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r= r, k = k)
      xxL <- lTableVec(xxMb, l = r*k)
      set.seed(1994)
      len <- length(xxL)
      for(indB in seq(1,B)){
        
        indsBst <- sample(seq(1, len), replace = T)
        bsts <- unlist(xxL[indsBst])# 
        bstSamp[indB] <- 
          tryCatch(
            gev_rl(Time, gev_mle_cpp_lvec(bsts)$par), 
            #gev_rl(Time, gev_mle_cpp(bsts)$par), 
            error = function(cond){
              cat("\n \n Erorr at indB = ", indB, "\n exact Error: ")
              message(conditionMessage(cond))
              return(c(NA,NA))}
          )
      }
    }else if(mthd == "db"){
      xxMb <- kMaxC(xx, r= r, k = 1)
      set.seed(1994)
      for(indB in seq(1,B)){
        bstSamp[indB] <-
          tryCatch(gev_rl(Time, evd::fgev(dbBootstrap(xxMb), std.err = F)$par), 
                   error = function(cond){
                     cat("\n \n Erorr at indB = ", indB, "\n exact Error: ")
                     message(conditionMessage(cond))
                     return(c(NA,NA))
                   })
      }
    }
    bstSamp <- na.omit(bstSamp)
    #print(quantile(bstSamp, c(0.01, 0.02, 0.03, 0.97, 0.98, 0.99)))
    quants <- quantile(bstSamp, c(1-niv/2, niv/2))
    
    
    #tDel <- difftime(Sys.time(), t0)
    #cat(format(tDe l))
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r= r, k = 0)
    }
    est <- gev_rl(Time, gev_mle_cpp(xxMb)$par)
    return(unname(c(est, 2*est - quants)))
    # return( list(chars = unname(c(est, 2*est - quants)), 
    #              variance = var(bstSamp)
    #                ))
  }
