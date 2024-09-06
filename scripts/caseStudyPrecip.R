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
# Aachen: downLink <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_00003_18910101_20110331_hist.zip"
downLink <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_01078_19520101_20231231_hist.zip"
downLink <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_01358_18900801_20231231_hist.zip"
downLink <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_02290_17810101_20231231_hist.zip"#Hohenpeißenberg
download.file(downLink, "data.zip")
# Aachen name: fileName <- "produkt_klima_tag_18910101_20110331_00003.txt" #Aachen
fileName <- "produkt_klima_tag_19520101_20231231_01078.txt" #Ddorf
fileName <- "produkt_klima_tag_18900801_20231231_01358.txt" # Fichterlberg
fileName <- "produkt_klima_tag_17810101_20231231_02290.txt" #Hohenpeißenberg
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
#data %>% filter(as.double(format(day, "%Y")) == 1927)
data <- data %>% filter(RSK >= 0) #starts at 1969
gapVec <- numeric(length(data$day))
for(ind in seq(1, gapVec %>% length())){
  #immer jahre reinpacken und dann table erstellen
  gapVec[ind] <- as.double(format(data$day[ind], "%Y"))
}
gapVec %>% table() #no observations from 1914 to 1930
rm(gapVec, ind)
#data <- data %>% filter(as.double(format(day, "%Y")) >=1930)
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
#data <- data %>% filter(as.double(format(day, "%Y")) >=1946) %>% 
data <- data %>%  mutate(year = as.double(format(day, "%Y")))
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
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "db",
    estim = DbArr[,1],
    lower = DbArr[,2],
    upper = DbArr[,3]
  ), 
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "cb",
    estim = CbArr[,1],
    lower = CbArr[,2],
    upper = CbArr[,3]
  ),
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
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
ciTibPlt$type <- 
  factor(ciTibPlt$type, levels = c("db", "cb"))
ciTibPlt <- 
  ciTibPlt %>% mutate(width = upper - lower, year = year - 40)
ciTibPlt$type <- 
  ciTibPlt$type %>% factor(
    levels = c("cb", "db"),
    labels = c("cb" = "cb(2)", "db" = "db")
  )
ciTibPlt
textSize <- 20
themePlot <- theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2),
                   strip.background = element_rect(color = "black", 
                                                   fill = "lightgrey", size = 0.2),
                   axis.title.x = element_text(size = textSize),
                   axis.title.y = element_text(size = textSize),
                   axis.text.y =element_text(size=textSize), 
                   axis.text.x =element_text(size=textSize),
                   strip.text.x = element_text(size = textSize),
                   strip.text.y = element_text(size = textSize),
                   plot.title = element_text(hjust = 0.5, size = textSize, 
                                             face = "bold"), 
                   #panel.background = element_rect(rgb(0.95, 0.95, 0.95, alpha = 1)),
                   legend.position = "right",
                   legend.title = element_text(size = textSize),
                   legend.text = element_text(size = textSize))
ownPalette <- #based on dark2
  c("cb(2)" = "#F8766D",  
    "cb(3)" = "#7CAE00",  
    "db" = "#00BFC4",  
    "sb" = "#C77CFF")


yearBounds <- 
  (ciTibPlt$year)[c(1, nY - winSize + 1)]
ciPlot <- 
  ciTibPlt %>% mutate(width = upper - lower) %>% 
  ggplot(aes(x = year))+
  geom_line(
    data = tibble(
      x = seq(yearBounds[1], yearBounds[2]),
      y = precDb[seq(1, yearBounds[2] - yearBounds[1] +1)]
    ),
    aes(x = x, y = (y - min(precDb))/5), col = "darkblue"
  )+
  geom_line(
    data = ciTibPlt,
    aes(x = year, y = estim), col = "black")+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4,
              linewidth = 0.1, col = "red")+
  facet_wrap(vars(type))+
  scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize + 1], 
                                  length.out = 4) %>% 
      round()
    )+
  scale_color_manual(
    values = ownPalette
  )+
  #facet_wrap(vars(ciType), 
  #           labeller = labeller(ciType = c(bstrDb = "db", bstr = "cb")))+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in",
    y = paste0("Rainfall heigth"),
    x = "Year"
  )

#new main:
twoYeaAvgs <- 
  array(dim = c(length(ciTibPlt$upper)/2, 2))
for(indY in seq(1, length(ciTibPlt$upper)/2)){
  twoYeaAvgs[indY,1] <- 
    ciTibPlt$lower[c(indY, indY+1)] %>% mean()
  twoYeaAvgs[indY,2] <- 
    ciTibPlt$upper[c(indY, indY+1)] %>% mean()
}

twoYeaAvgs <-

ciPlot <- 
  ciTibPlt %>% 
  ggplot(aes(x = year))+
  geom_line(
    data = tibble(
      x = seq(yearBounds[1], yearBounds[2])+40,
      y = precDb[seq(1, yearBounds[2] - yearBounds[1] +1)]
    ),
    aes(x = x, y = (y - min(precDb))/5), col = "black"
  )+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.4,
              linewidth = 0.1, col = "red")+
  geom_line(
    data = ciTibPlt,
    aes(x = year, y = estim, col = type), linewidth = 1)+
  scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize + 1], 
          length.out = 4) %>% 
      round()
  )+
  scale_color_manual(
    values = ownPalette
  )+
  #facet_wrap(vars(ciType), 
  #           labeller = labeller(ciType = c(bstrDb = "db", bstr = "cb")))+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in",
    y = paste0("Precipitation"),
    x = "Year",
    col = "Bootstrap:",
    fill = "Bootstrap:"
  )
ciPlot


#explore begin
ciTibPlt %>% mutate(width = upper - lower, year = year + 40) %>% 
  ggplot(aes(x = year, y = width, col = type))+
  geom_line()+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize + 1], 
          length.out = 4) %>% 
      round()) + 40
  )+
  #facet_wrap(vars(ciType), 
  #           labeller = labeller(ciType = c(bstrDb = "db", bstr = "cb")))+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("CI Width"),
    x = "Year"
  )


#explore end


if(F){
  ggsave(plot = ciPlot, filename = here("results/plotCaseStudyCbandsMain.pdf"),
         device = "pdf", width = 10, height = 4)
}

#DEPRECATED
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
ciTib <- ciTib %>% mutate(
  width = upper - lower) %>% group_by(type) %>% summarise(avgWidth = median(width))

#plotting of confidence intervals

ciTibPlt <- ciTib
ciTibPlt <- ciTibPlt %>% filter(type != "sb")

ciTibPlt %>% ggplot(aes(x = year))+
  geom_line(aes(x = year, y = estim), col = "black")+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, linewidth = 0.4)+
  facet_wrap(vars(type))

