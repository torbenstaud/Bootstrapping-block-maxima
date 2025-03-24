#load coeff_rl vector

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
  
  # Versuch mit "BFGS"
  resultBfgs <- try({
    optim(
      par = start,
      fn = function(theta) neg_log_likelihood_gev_lvec(theta, data),
      method = "BFGS",
      #lower = c(-10, 0.01, -0.9),
      #upper = c(50, 100, 0.9),
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


ciCircmax <- 
  function(xx, Time = 100, r, k, niv = 0.05, B = 10^3,
           mthd = "cb"){
    #t0 <- Sys.time()
    array_bst <- array(dim = c(B+1,3))
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r= r, k = k)
      xxL <- lTableVec(xxMb, l = r*k)
      set.seed(1994)
      len <- length(xxL)
      for(indB in seq(1,B)){
        
        indsBst <- sample(seq(1, len), replace = T)
        bsts <- unlist(xxL[indsBst])# 
        array_bst[1+indB,] <- 
          tryCatch(
            gev_mle_cpp_lvec(bsts)$par, 
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
        array_bst[1+indB,] <-
          tryCatch(evd::fgev(dbBootstrap(xxMb), std.err = F)$par, 
                   error = function(cond){
                     cat("\n \n Erorr at indB = ", indB, "\n exact Error: ")
                     message(conditionMessage(cond))
                     return(c(NA,NA))
                   })
      }
    }
    
    
    #tDel <- difftime(Sys.time(), t0)
    #cat(format(tDe l))
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r= r, k = 0)
    }
    array_bst[1,] <- gev_mle_cpp(xxMb)$par
    
    return(array_bst)
    # return( list(chars = unname(c(est, 2*est - quants)), 
    #              variance = var(bstSamp)
    #                ))
  }


ciCircmaxRl <- 
  function(xx, Time = 100, r, k, niv = 0.05, B = 10^3, topcutoffGam = 0.5, 
           botcutoffGam = -0.3, onesided = F,
           mthd = "cb"){
    #t0 <- Sys.time()
    m <- as.integer(length(xx)/r)
    bstSamp <- rep(NA, B)
    theta_bst <- c(100,100,100)
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r= r, k = k)
      xxL <- lTableVec(xxMb, l = r*k)
      set.seed(1994)
      len <- length(xxL)
      
      for(indB in seq(1,B)){
        while(!(theta_bst[3] <= topcutoffGam & theta_bst[3] >= botcutoffGam)){
          indsBst <- sample(seq(1, len), replace = T)
          bsts <- unlist(xxL[indsBst])# 
          theta_bst <- gev_mle_cpp_lvec(bsts)$par
        }
        
        bstSamp[indB] <-  gev_rl(Time,theta_bst)
        theta_bst[3] <- 100
        # if(theta[3] <= topcutoffGam & theta[3] >= botcutoffGam){
        #   
        #     }
      }
    }else if(mthd == "db"){
      xxMb <- kMaxC(xx, r= r, k = 1)
      set.seed(1994)
      for(indB in seq(1,B)){
        while(!(theta_bst[3] <= topcutoffGam & theta_bst[3] >= botcutoffGam)){
          indBst <- sample(seq(1, m), replace = T)
          bsts <- xxMb[indBst]
          theta_bst <- evd::fgev(bsts, std.err = F)$par
          
        }
        bstSamp[indB] <-  gev_rl(Time,theta_bst)
        # if(theta[3] <= topcutoffGam & theta[3] >= botcutoffGam){
        #   
        # }
        theta_bst[3] <- 100
      }
    }
    # pct_omit <- sum(is.na(bstSamp))*100/B
    # bstSamp <- na.omit(bstSamp) 
    # cat(sprintf("Es wurden %.2f Prozent der Bootstraps verworfen im Verfahren %s.\n", 
    #             pct_omit, mthd))
    #print(quantile(bstSamp, c(0.01, 0.02, 0.03, 0.97, 0.98, 0.99)))
    quants <- quantile(bstSamp, c(1-niv/2, niv/2), na.rm = T)
    
    #tDel <- difftime(Sys.time(), t0)
    #cat(format(tDe l))
    theta_bst <- gev_mle_cpp(xxMb)$par
    est_bst <- gev_rl(Time, theta_bst)
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r=r, k = 0)
    }
    theta_est <- gev_mle_cpp(xxMb)$par
    
    est <- gev_rl(Time, theta_est)
    # ##hier noch factor correction hin: basierend auf m und \hat gamma
    # #quants_corr <- quants*c_corr(m, theta_est[3]) #c_corr noch defn
    # c <- coeff_rl[1] + theta_est[3]*coeff_rl[2] + length(xx)/r*coeff_rl[3]
    # quants_corr <- (est_bst - quants)*c
    #cat( sprintf("CI ist: [%s].\n",est + quants_corr))
    return(unname(c(est, est_bst + est - quants)))
    # return( list(chars = unname(c(est, 2*est - quants)), 
    #              variance = var(bstSamp)
    #                ))
  }

ciCircmaxMean <- 
  function(xx, Time = 100, r, k, niv = 0.05, B = 10^3, topcutoffGam = 100.5, 
           botcutoffGam = -100.5, mthd = "cb", onesided = F){
    #t0 <- Sys.time()
    m <- as.integer(length(xx)/r)
    bstSamp <- rep(NA, B)
    theta_bst <- c(100,100,100)
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r= r, k = k)
      xxL <- lTableVec(xxMb, l = r*k)
      set.seed(1994)
      len <- length(xxL)
      
      for(indB in seq(1,B)){
        counter <- 0
        indsBst <- sample(seq(1, len), replace = T)
        bsts <- unlist(xxL[indsBst])#
        # while(!(theta_bst[3] <= topcutoffGam & theta_bst[3] >= botcutoffGam)){
        #   indsBst <- sample(seq(1, len), replace = T)
        #   bsts <- unlist(xxL[indsBst])# 
        #   theta_bst <- gev_mle_cpp_lvec(bsts)$par
        #   #counter <- counter +1
        #   #cat(counter)
        # }
        
        bstSamp[indB] <-  meanCTabVec(bsts)
        theta_bst[3] <- 100
        # if(theta[3] <= topcutoffGam & theta[3] >= botcutoffGam){
        #   
        #     }
      }
    }else if(mthd == "db"){
      xxMb <- kMaxC(xx, r= r, k = 1)
      set.seed(1994)
      for(indB in seq(1,B)){
        indBst <- sample(seq(1, m), replace = T)
        bsts <- xxMb[indBst]
        # while(!(theta_bst[3] <= topcutoffGam & theta_bst[3] >= botcutoffGam)){
        #   indBst <- sample(seq(1, m), replace = T)
        #   bsts <- xxMb[indBst]
        #   theta_bst <- evd::fgev(bsts, std.err = F)$par
        #   
        # }
        bstSamp[indB] <-  mean(bsts)
        theta_bst[3] <- 100
      }
    }
    # pct_omit <- sum(is.na(bstSamp))*100/B
    # bstSamp <- na.omit(bstSamp) 
    # cat(sprintf("Es wurden %.2f Prozent der Bootstraps verworfen im Verfahren %s.\n", 
    #             pct_omit, mthd))
    #print(quantile(bstSamp, c(0.01, 0.02, 0.03, 0.97, 0.98, 0.99)))
    if(onesided){
      quants <- c(-Inf, quantile(bstSamp, niv/2, na.rm = T))
    }else{
      quants <- quantile(bstSamp, c(1-niv/2, niv/2), na.rm = T)
    }
    
    #tDel <- difftime(Sys.time(), t0)
    #cat(format(tDe l))
    est_bst <- mean(xxMb)
    if(mthd == "cb"){
      xxMb <- kMaxC(xx, r=r, k = 0)
    }
    theta_est <- gev_mle_cpp(xxMb)$par
    est <- mean(xxMb)
    # ##hier noch factor correction hin: basierend auf m und \hat gamma
    # #quants_corr <- quants*c_corr(m, theta_est[3]) #c_corr noch defn
    # c <- coeff_mean[1] + theta_est[3]*coeff_mean[2] + length(xx)/r*coeff_mean[3]
    # quants_corr <- (est_bst - quants)*c
    #cat( sprintf("CI ist: [%s].\n",est + quants_corr))
    return(unname(c(est, est_bst + est - quants)))
    #return(unname(c(est, est_bst + est - quants,bstSamp))) #debug
    
    
    # return( list(chars = unname(c(est, 2*est - quants)), 
    #              variance = var(bstSamp)
    #                ))
  }
