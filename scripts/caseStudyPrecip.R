#load libraries----
library(readr) #read and download files
library(tidyr) #data wrangling
library(dplyr) #data wrangling
library(ggplot2) #plotting
library(here) #relative paths
library(Rcpp) #cpp api
#source files----
source(here("src/0source.R"))
sourceCpp(here("src/0sourceC.cpp"))
source(here("src/sourceCase.R"))
source(here("src/theme.R"))


#download precipitation data (location: Aachen, Germany) and read it
##link leads to an zip-archive
downLink <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_00003_18910101_20110331_hist.zip"
download.file(downLink, "data.zip")
fileName <- "produkt_klima_tag_18910101_20110331_00003.txt"
dataRaw <- unz("data.zip", fileName) %>% read_delim(delim=";", trim_ws = T)
##delete the downloaded file
file.remove("data.zip")
if(!file.exists("data.zip")){
  print("File succesfully deleted!")
}else{
  print("File could not be deleted!")
}
rm(downLink, fileName)
##only need precipitation height in mm (RSK) and date of measuring (MESS_DATUM)
data <- dataRaw %>% select(MESS_DATUM, RSK) %>% 
  rename(day = "MESS_DATUM") %>% 
  mutate(day = as.Date(as.character(day), format = "%Y%m%d"))
#rm(dataRaw)
data
## missing or absurd data?
data %>% na.omit()
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
data %>% filter(as.double(format(day, "%Y")) == 1927)
gapVec <- numeric(length(data$day))
for(ind in seq(1, gapVec %>% length())){
  #immer jahre reinpacken und dann table erstellen
  gapVec[ind] <- as.double(format(data$day[ind], "%Y"))
}
gapVec %>% table() #no observations from 1914 to 1930
rm(gapVec, ind)
data <- data %>% filter(as.double(format(day, "%Y")) >=1930)
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
data %>% filter(RSK < 0) #false data in 1945
data <- data %>% filter(as.double(format(day, "%Y")) >=1946) %>% 
  mutate(year = as.double(format(day, "%Y")))
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
# data <- data %>% bind_rows(
#   tibble(day = data$day + length(data$day), RSK = data$RSK) %>% 
#     mutate(year = as.double(as.double(format(day, "%Y"))))
# )
# data <- data %>% bind_rows(
#   tibble(day = data$day + length(data$day), RSK = data$RSK) %>% 
#     mutate(year = as.double(as.double(format(day, "%Y"))))
# )
nY <- max(data$year) - min(data$year) + 1
## -> left with 66 full years of daily observations

## diagnostics: are disjoint block maxima over a year approx independent?
dbmsFull <- data %>% group_by(year) %>% summarize(blockmax = max(RSK))
dbmsFull$blockmax %>% acf() #no significant evidence against independence
#plan: rolling window over 40 years and 1 year jumps
winSize <-  40
blockSize <- 365 #one year blocks

##vector of disjoint yearly precipitation maxima
precDb <- numeric(nY)
yearVec <- (data$year) %>% unique()
for(indY in seq_along(precDb)){
  precDb[indY] <- 
    (data %>% filter(year == yearVec[indY]))$RSK %>% max()
}
##calculate estimators and bootstrap CIs
DbArr <- array( dim = c(nY - winSize + 1, 3))
SbArr <- array( dim = c(nY - winSize + 1, 3))
CbArr <- array( dim = c(nY - winSize + 1, 3))
bstVars <- array(dim = c(nY - winSize + 1, 2))
for(indW in seq(1, nY - winSize + 1)){
  dataWin <- (data %>% filter(
    year %in% seq(yearVec[indW], yearVec[indW] + 39)
  ))$RSK
  DbArr[indW,1] <- dataWin %>% 
    kMaxC(r = blockSize, k = 1) %>% 
    mean() #for mean study
  ciDbDat <- 
    ciCircmaxMean(dataWin, B = 1000, r = blockSize, k = 1, mthd = "db")
  DbArr[indW,c (2,3)] <- ciDbDat[[1]]
  bstVars[indW,1] <- ciDbDat[[2]]
  SbArr[indW,1] <- (dataWin) %>% 
    kMaxC(r = blockSize, k = 0) %>% 
    mean()
    
  CbArr[indW,1] <- (dataWin) %>% 
    kMaxC(r = blockSize, k = 2) %>% 
    mean()
  ciCbDat <- 
    ciCircmaxMean(dataWin, B = 1000, r = blockSize, k = 2, mthd = "cb")
  CbArr[indW,c (2,3)] <- ciCbDat[[1]]
  bstVars[indW,2] <- ciCbDat[[2]]
}
# create tibble with data
ciTib <- bind_rows(
  tibble(
    year = seq(1,nY - winSize + 1),
    type = "db",
    estim = DbArr[,1],
    lower = DbArr[,2],
    upper = DbArr[,3]
  ), 
  tibble(
    year = seq(1,nY - winSize + 1),
    type = "cb",
    estim = CbArr[,1],
    lower = CbArr[,2],
    upper = CbArr[,3]
  ),
  tibble(
    year = seq(1,nY - winSize + 1),
    type = "sb",
    estim = SbArr[,1],
    lower = SbArr[,2],
    upper = SbArr[,3]
  )
  
)
ciTib %>% mutate(
  width = upper - lower
) %>% View()

avgWidths <- (ciTib %>% mutate(
  width = upper - lower) %>% group_by(type) %>% 
  summarise(avgWidth = median(width)) %>% 
  filter(type != "sb") %>% select(avgWidth))[[1]]
avgWidths[2]/avgWidths[1]
#plotting of confidence intervals

ciTibPlt <- ciTib
ciTibPlt <- ciTibPlt %>% filter(type != "sb")

ciTibPlt %>% ggplot(aes(x = year))+
  geom_line(aes(x = year, y = estim), col = "black")+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, linewidth = 0.4)+
  facet_wrap(vars(type))

#analysis for the estimation of Prob(M_r <= threshold)
##threshold for a yearly Maximum will be set as the empirical median of the DB maxima
thresh <- precDb %>% quantile(0.5) %>% unname()
threshInd <- function(xx){
  ifelse(xx <= thresh, 1, 0)
}
DbArr <- array( dim = c(nY - winSize + 1, 3))
SbArr <- array( dim = c(nY - winSize + 1, 3))
CbArr <- array( dim = c(nY - winSize + 1, 3))
bstVars <- array(dim = c(nY - winSize + 1, 2))
for(indW in seq(1, nY - winSize + 1)){
  dataWin <- (data %>% filter(
    year %in% seq(yearVec[indW], yearVec[indW] + 39)
  ))$RSK
  DbArr[indW,1] <- dataWin %>% 
    kMaxC(r = blockSize, k = 1) %>% 
    threshInd() %>% mean()
  #mean() #for mean study
  ciDbDat <- 
    ciCircmaxProb(dataWin, B = 1000, r = blockSize, k = 1, mthd = "db", 
                  tresh = thresh)
  #ciCircmaxMean(dataWin, B = 250, r = blockSize, k = 1, mthd = "db")
  DbArr[indW,c (2,3)] <- ciDbDat[[1]]
  bstVars[indW,1] <- ciDbDat[[2]]
  SbArr[indW,1] <- (dataWin) %>% 
    kMaxC(r = blockSize, k = 0) %>% 
    threshInd() %>% mean()
  #var()
  CbArr[indW,1] <- (dataWin) %>% 
    kMaxC(r = blockSize, k = 2) %>% 
    threshInd() %>% mean()
  #var()
  ciCbDat <- ciCircmaxProb(dataWin, B = 1000, r = blockSize, k = 2, mthd = "cb", 
                           tresh = thresh)
  #ciCircmaxMean(dataWin, B = 250, r = blockSize, k = 2, mthd = "cb")
  CbArr[indW,c (2,3)] <- ciCbDat[[1]]
  bstVars[indW,2] <- ciCbDat[[2]]
}
DbArr
SbArr
CbArr
# create tibble with data
ciTib <- bind_rows(
  tibble(
    year = seq(1,nY - winSize + 1),
    type = "db",
    estim = DbArr[,1],
    lower = DbArr[,2],
    upper = DbArr[,3]
  ), 
  tibble(
    year = seq(1,nY - winSize + 1),
    type = "cb",
    estim = CbArr[,1],
    lower = CbArr[,2],
    upper = CbArr[,3]
  ),
  tibble(
    year = seq(1,nY - winSize + 1),
    type = "sb",
    estim = SbArr[,1],
    lower = SbArr[,2],
    upper = SbArr[,3]
  )
  
)
ciTib %>% mutate(
  width = upper - lower
) %>% View()
ciTib %>% mutate(
  width = upper - lower) %>% group_by(type) %>% summarise(avgWidth = median(width))
#plotting of confidence intervals

ciTibPlt <- ciTib
ciTibPlt <- ciTibPlt %>% filter(type != "sb")

ciTibPlt %>% ggplot(aes(x = year))+
  geom_line(aes(x = year, y = estim), col = "black")+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, linewidth = 0.4)+
  facet_wrap(vars(type))

