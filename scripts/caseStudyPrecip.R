#load libraries
library(readr) #read and download files
library(tidyr) #data wrangling
library(dplyr) #data wrangling
library(ggplot2) #plotting

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
data <- data %>% filter(as.double(format(day, "%Y")) >=1946)
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
nY <- 66
## -> left with 66 full years of daily observations

## diagnostics: are disjoint block maxima over a year approx independent?
dbmsFull <- data %>% group_by(year) %>% summarize(blockmax = max(RSK))
dbmsFull$blockmax %>% acf() #no significant evidence against independence
#plan: rolling window over 40 years and 1 year jumps
winSize <- 40
data <- data %>% mutate(year = as.double(format(day, "%Y")))

nY
##vector of disjoint yearly precipitation maxima
precDb <- numeric(nY)
yearVec <- (data$year) %>% unique()
for(indY in seq_along(precDb)){
  precDb[indY] <- 
    (data %>% filter(year == yearVec[indY]))$RSK %>% max()
}
##threshold for a yearly Maximum will be set as the empirical median of the DB maxima
thresh <- precDb %>% quantile(0.5) %>% unname()
threshInd <- function(xx){
  ifelse(xx <= thresh, 1, 0)
}
estDb <- numeric(nY - winSize + 1)
estSb <- numeric(nY - winSize + 1)
estCb <- numeric(nY - winSize +1)

for(indW in seq_along(estDb)){
  dataWin <- data %>% filter(
    year %in% seq(yearVec[indW], yearVec[indW] + 39)
  )
  estDb[indW] <- (dataWin$RSK) %>% 
    RcppRoll::roll_max(n = 365, by = 365, na.rm = T) %>% 
    threshInd() %>% mean()
  estSb[indW] <- (dataWin$RSK) %>% 
    RcppRoll::roll_max(n = 365, na.rm = T) %>% 
    threshInd() %>% mean()
  estCb[indW] <- (dataWin$RSK) %>% 
    RcppRoll::roll_max(n = 365, na.rm = T) %>% 
    threshInd() %>% mean()
}
estDb
estSb
