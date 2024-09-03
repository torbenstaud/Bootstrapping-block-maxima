library(evd)
#vergleichbare n,r wie in case study
N <- 160
n <- 2500
r <- 62
niv <- 0.05
B <- 10^3
#treshold durch median eines dbMax samples tarieren
tresh <- rgev(n,loc = 1.665e-2, scale = 5.553e-3, shape = 1.356e-6) %>% 
  kMaxC(r, 1) %>% quantile(0.5)
## initialize data structures
estMat <- array(dim = c(N,4))
cbWidthsMat <- array(dim = c(N,3))
ciArray <- array(dim = c(N, 3 , 2))
for(indN in seq(1,N)){
  set.seed(indN)
  xxs <- rgev(n,loc = 1.665e-2, scale = 5.553e-3, shape = 1.356e-6)
  treshInd <- function(xx){
    ifelse(xx <= tresh, 1, 0)
  }
  estMat[indN, 1] <- kMaxC(xxs, r, 0) %>% treshInd() %>% mean()
  estMat[indN, 2] <- kMaxC(xxs, r, 1) %>% treshInd() %>% mean()
  estMat[indN, 3] <- kMaxC(xxs, r, 2) %>% treshInd() %>% mean()
  estMat[indN, 4] <- kMaxC(xxs, r, 3) %>% treshInd() %>% mean()
  #ci widths für db, cb(2), cb(3) berechnen basierend auf bootstrap
  ciTmp <- ciCircmaxProb(xxs, r, 1, B = B, tresh = tresh, mthd = "db")[[1]]
  cbWidthsMat[indN, 1] <- 
    ciTmp %>% diff()
  ciArray[indN, 1, ] <- 2*estMat[indN, 2] - ciTmp #db est - db quantiles
  
  ciTmp <- ciCircmaxProb(xxs, r, 2, B = B, tresh = tresh)[[1]] 
  cbWidthsMat[indN, 2] <- 
    ciTmp %>% diff()
  ciArray[indN, 2,] <- 2*estMat[indN, 1] - ciTmp #sb est - cb(2) quantiles
  
  ciTmp <- ciCircmaxProb(xxs, r, 3, B = B, tresh = tresh)[[1]]  
  cbWidthsMat[indN, 3] <- 
    ciTmp %>% diff()
  ciArray[indN, 3, ] <- 2*estMat[indN, 1] - ciTmp #sb est - cb(3) quantiles
}
estMat %>% apply(c(2), mean) #große Unterschiede? (sollten keine sein)
vars <- estMat %>% apply(c(2), var) #empVar der schätzer (sollte sb < cb(3) < cb(2) < db())
vars[2]/vars #ratios von drüber
(cbWidthsMat[,1]/cbWidthsMat[,2]) %>% mean() #>1 heißt: cb(2) klappt besser
##visuell mit ggplot
plotTib <- 
  bind_rows(
    tibble(day = seq(1,N),
           type = "db",
           lower = ciArray[,1,1],
           upper = ciArray[,1,2],
           est = estMat[,1]),
    tibble(day = seq(1,N),
           type = "cb(2)",
           lower = ciArray[,2,1],
           upper = ciArray[,2,2],
           est = estMat[,2]),
    tibble(day = seq(1,N),
           type = "cb(3)",
           lower = ciArray[,3, 1],
           upper = ciArray[,3, 2],
           est = estMat[,2])
  )
plotTib$type <- plotTib$type %>% factor(
  levels = c("db", "cb(2)", "cb(3)"),
  labels = c("db", "cb", "cb(3)")
)
plotTib %>% filter(type != "cb(3)", day %%2 == 0) %>% 
  ggplot(aes(x=day, y = lower))+
  geom_ribbon(
    aes(ymin = lower, ymax = upper), col = "red", alpha = 0.4, linewidth = 0.1
  )+
  geom_line(aes(x= day, y = est), col = "black")+
  facet_wrap(vars(type))+
  scale_y_continuous(limits = c(0,1))

  
#now the study with a sliding window as in the case study
N <- 160
n <- 2500
r <- 62
niv <- 0.05
B <- 1000

## initialize data structures
estMat <- array(dim = c(N,4))
cbWidthsMat <- array(dim = c(N,3))
ciArray <- array(dim = c(N, 3 , 2))
##create data (die weirden parameter kommen von einem ML GEV fit der case study daten)
set.seed(1234)
data <- 
  rgev(r*4*50,loc = 1.665e-2, scale = 5.553e-3, shape = 1.356e-6)

## median des GESAMTEN db samples als tresholds
tresh <- data %>% kMaxC(r,1) %>% quantile(0.5)
##die 50Jahre je Quartal durchsliden und immer ein 10 Jahresfenster pro indN angucken
for(indN in seq(1,N)){ 
  xxs <- #data for 10 years window sliding by quarters
    data[seq( (indN - 1)*r+1, (indN +39)*62 )]
  treshInd <- function(xx){
    ifelse(xx <= tresh, 1, 0)
  }
  #Schätzer (slid, db, cb(2), cb(3)) auswerten
  estMat[indN, 1] <- kMaxC(xxs, r, 0) %>% treshInd() %>% mean()
  estMat[indN, 2] <- kMaxC(xxs, r, 1) %>% treshInd() %>% mean()
  estMat[indN, 3] <- kMaxC(xxs, r, 2) %>% treshInd() %>% mean()
  estMat[indN, 4] <- kMaxC(xxs, r, 3) %>% treshInd() %>% mean()
  
  #ci widths für db, cb(2), cb(3) berechnen basierend auf bootstrap
  ciTmp <- ciCircmaxProb(xxs, r, 1, B = B, tresh = tresh, mthd = "db")[[1]]
  cbWidthsMat[indN, 1] <- 
    ciTmp %>% diff()
  ciArray[indN, 1, ] <- ciTmp
  
  ciTmp <- ciCircmaxProb(xxs, r, 2, B = B, tresh = tresh)[[1]] 
  cbWidthsMat[indN, 2] <- 
    ciTmp %>% diff()
  ciArray[indN, 2,] <- ciTmp
  
  ciTmp <- ciCircmaxProb(xxs, r, 3, B = B, tresh = tresh)[[1]]  
  cbWidthsMat[indN, 3] <- 
    ciTmp %>% diff()
  ciArray[indN, 3, ] <- ciTmp
}  
estMat %>% apply(c(2), mean) #große Unterschiede? (sollten keine sein)
vars <- estMat %>% apply(c(2), var) #empVar der schätzer (sollte sb < cb(3) < cb(2) < db())
vars[2]/vars #ratios von drüber
(cbWidthsMat[,1]/cbWidthsMat[,2]) %>% mean() #>1 heißt: cb(2) klappt besser
#(cbWidthsMat[,1]/cbWidthsMat[,2]) %>% ts.plot()

##visuell mit ggplot
plotTib <- 
  bind_rows(
  tibble(day = seq(1,N),
         type = "db",
         lower = ciArray[,1,1],
         upper = ciArray[,1,2],
         est = estMat[,2]),
  tibble(day = seq(1,N),
         type = "cb(2)",
         lower = ciArray[,2,1],
         upper = ciArray[,2,2],
         est = estMat[,3]),
  tibble(day = seq(1,N),
         type = "cb(3)",
         lower = ciArray[,3, 1],
         upper = ciArray[,3, 2],
         est = estMat[,4])
  )
plotTib$type <- plotTib$type %>% factor(
  levels = c("db", "cb(2)", "cb(3)"),
  labels = c("db", "cb", "cb(3)")
)
plotTib %>% filter(type != "cb(2)124") %>% 
  ggplot(aes(x=day, y = lower))+
  geom_line(aes(x=day, y = est), col = "black")+
  geom_ribbon(
    aes(ymin = lower, ymax = upper), col = "red", alpha = 0.4, linewidth = 0.1
  )+
  facet_wrap(vars(type))+
  scale_y_continuous(limits = c(0,1))

ts.plot(cbWidthsMat[,1]) #db
ts.plot(cbWidthsMat[,2]) #cb(2)
ts.plot(cbWidthsMat[,3]) #cb(3)


#they bands are way to irregular: what happens with iid normal bands for mean estim
N <- 160
n <- 200 #ceiling(2500/60) #approx 42
estVec <- numeric(N)
ciArr <- array(dim = c(N,2))
normBstr <- function(xxs, B = 10^3){
  est <- mean(xxs)
  res <- numeric(2)
  bstr <- numeric(B)
  for(indB in seq(1,B)){
    set.seed(indB)
    bstr[indB] <- 
      mean(dbBootstrap(xxs))
  }
  res <- 2*est - quantile(bstr, c(1-0.05/2, 0.05/2))
}
for(indN in seq(1,N)){
  set.seed(indN)
  xxs <- rnorm(n)
  estVec[indN] <- mean(xxs)
  ciArr[indN,] <- normBstr(xxs)
}
plotTib <- 
  tibble(day = 1:160,
         est = estVec,
         lower = ciArr[,1],
         upper = ciArr[,2])
plotTib %>% ggplot(aes(x = day ,y = est)) +
  geom_line(col = "black"
  )+
  geom_ribbon(
    aes(ymin = lower, ymax = upper), alpha = 0.4
  )

#clean up  
rm(testSamp)