#load libraries----
library(readr) #read and download files
library(tidyr) #data wrangling
library(dplyr) #data wrangling
library(ggplot2) #plotting
library(ggpubr)#arranging 
library(here) #relative paths
library(Rcpp) #cpp api
#source files----
source(here("src/0source.R"))
source(here("src/0source_ts.R"))
sourceCpp(here("src/0sourceC.cpp"))
sourceCpp(here("src/0sourceCn.cpp"))
source(here("src/sourceCase.R"))
source(here("src/theme.R"))


#download precipitation data (location: Aachen, Germany) and read it
##link leads to an zip-archive
downLink <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_02290_17810101_20231231_hist.zip"#Hohenpeißenberg
download.file(downLink, "data.zip")
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
data <- dataRaw %>% select(MESS_DATUM, RSK, TXK) %>% 
  rename(day = "MESS_DATUM") %>% 
  mutate(day = as.Date(as.character(day), format = "%Y%m%d"))
#rm(dataRaw)
data
## missing or absurd data?
data %>% na.omit()
data %>% ggplot(aes(x = day, y = TXK))+
  geom_line()
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
#data %>% filter(as.double(format(day, "%Y")) == 1927)
data <- data %>% filter(RSK >= 0, TXK > -999) #starts at 1969
gapVec <- numeric(length(data$day))
for(ind in seq(1, gapVec %>% length())){
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
data %>% ggplot(aes(x = day, y = TXK))+
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
dbmsFull <- data %>% group_by(year) %>% summarize(blockmax = max(TXK))
dbmsFull$blockmax %>% acf() #no significant evidence against independence
#plan: rolling window over 40 years and 1 year jumps
winSize <-  40
blockSize <- 365 #one year blocks

##vector of disjoint yearly precipitation/temp maxima
precDb <- numeric(nY)
tempDb <- numeric(nY)
yearVec <- (data$year) %>% unique()
for(indY in seq_along(precDb)){
  precDb[indY] <- 
    (data %>% filter(year == yearVec[indY]))$RSK %>% max()
  tempDb[indY] <- 
    (data %>% filter(year == yearVec[indY]))$TXK %>% max()
}
##calculate estimators and bootstrap CIs
DbArr <- array( dim = c(nY - winSize + 1, 3))
SbArr <- array( dim = c(nY - winSize + 1, 3))
CbArr <- array( dim = c(nY - winSize + 1, 3))
bstVars <- array(dim = c(nY - winSize + 1, 2))

if(F){#17mins with B = 10**3
  B <-  10**3
  t0 <- Sys.time()
  for(indW in seq(1, nY - winSize + 1)){
    dataWin <- (data %>% filter(
      year %in% seq(yearVec[indW], yearVec[indW] + 39)
    ))$RSK
    ciDbDat <- 
      ciCircmaxRl(dataWin, B = B, r = blockSize, k = 1, mthd = "db", onesided = F,
                  botcutoffGam = -0.3, topcutoffGam = 0.5) #werte für die cutoffs, da die shape über den ganzen Zeitraum als ~0.1 geschätzt wird und wir ein 0.4 band drum legen
    DbArr[indW,c(1, 2,3)] <- ciDbDat #est, lower and upper
    
    if(max(abs(ciDbDat)) >= 100 ){
      sprintf("Error at indW = %g für db", indW)
    }
    
    # ciSbDat <- 
    #   ciCircmaxRl(dataWin, B = B, r = blockSize, k = 1, mthd = "cb")
    # SbArr[indW,c(1, 2,3)] <- ciSbDat[[1]]
    # bstVars[indW,2] <- ciCbDat[[2]]
    # print("sb")#debug
    
    ciCbDat <- 
      ciCircmaxRl(dataWin, B = B, r = blockSize, k = 2, mthd = "cb", onesided = F)
    CbArr[indW,c(1, 2,3)] <- ciCbDat #est, lower and upper
    
    if(max(abs(ciCbDat)) >= 100 ){
      sprintf("Error at indW = %g für cb", indW)
    }
    
    if(indW %% 10 == 0) {print(indW)}
  }
  tDel <- difftime(Sys.time(), t0)
  cat(format(tDel))
  save(DbArr, CbArr, file = file.path("data", "array_rl_cs"))
}
load(file = file.path("data", "array_rl_cs_-0.30.5")) # DbArr, CbArr
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
    estim = CbArr[,1],
    lower = SbArr[,2],
    upper = SbArr[,3]
  )
) %>% mutate(
  year = year + 40
)
ciTib %>% mutate(
  width = upper - lower
) 

avgWidths <- (ciTib %>% mutate(
  width = upper - lower) %>% group_by(type) %>% 
    summarise(avgWidth = median(width)) %>% 
    filter(type != "sb") %>% select(avgWidth))[[1]]
avgWidths[2]/avgWidths[1]
##averaging of two year window
twoYeaAvgs <- 
  array(dim = c(length(ciTib$upper)/2, 2))
tmpTib <- tibble()
for(indT in seq(1,3)){
  indTT <- (ciTib$type %>% unique())[indT]
  updTib <- 
    (ciTib %>% filter(type == indTT))[,c(3,4,5)] %>% 
    apply(c(2), function(xx){stats::filter(xx, filter = c(1/2,1/2), sides = 1)}) %>% 
    as_tibble()
  updTib <- bind_cols(
    (ciTib %>% filter(type == indTT))[, c(1,2)],
    updTib
  )
  tmpTib <- tmpTib %>% bind_rows(updTib)
}
avgCiTib <- tmpTib %>% filter(year != 1919)
rm(tmpTib)



#plotting of confidence intervals

ciTibPlt <- avgCiTib
##cb estim substituted by sb estimator
ciTibPlt <- ciTibPlt %>%  bind_cols(
  sbEstim = rep((ciTibPlt %>% filter(type == "sb"))$estim, 3)
) %>% mutate(
  estim = ifelse(type == "cb", sbEstim, estim)
) %>% select(-c(sbEstim))

ciTibPlt <- ciTibPlt %>% filter(type != "sb")
ciTibPlt$type <- 
  factor(ciTibPlt$type, levels = c("db", "cb"))
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
  (ciTibPlt$year)[c(1, nY - winSize)]


#new main:

##averaged db block maxima
avgPrecDb <- (precDb %>% 
                stats::filter(sides = 1, filter = c(1/2,1/2)))[-1]

ciPlot <- 
  ciTibPlt %>%  
  ggplot(aes(x = year))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.4,
              linewidth = 0.1, col = "black")+
  geom_line(
    data = ciTibPlt,
    aes(x = year, y = estim, col = type), linewidth = 1)+
 # scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
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
ciTibPlt <- ciTibPlt %>% mutate(width = upper - lower) 
ciWidthPlt <- 
  ciTibPlt %>% 
  ggplot(aes(x = year, y = width, col = type))+
  geom_line(linewidth = 1.1)+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
           length.out = 4) %>% 
         round()) 
  )+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("CI width"),
    x = "Year"
  )
ciWidthPlt


fullCiPlot <- 
  ggarrange(
    ciPlot, ciWidthPlt, common.legend = T, legend = "right"
  )
fullCiPlot

#explore end
if(F){
  ggsave(plot = fullCiPlot, filename = here("results/plotCaseStudyCbandsMain.pdf"),
         device = "pdf", width = 10, height = 5) #ursprl 4
}
if(F){
  ggsave(plot = fullCiPlot2, filename = here("results/plotCaseStudyCbandsRel.pdf"),
         device = "pdf", width = 10, height = 5) #ursprl 4
}

#Jrssb

textSize <- 15
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

ciPlotJrssb <- 
  ciTibPlt %>% 
  ggplot(aes(x = year))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type, col = type), alpha = 0.4,
              linewidth = 0.5)+
  geom_line(
    data = ciTibPlt,
    aes(x = year, y = estim, col = type), linewidth = 1.1)+
  #scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
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
ciPlotJrssb

ciWidthPltJrssb <- 
  ciTibPlt %>% mutate(width = upper - lower) %>%  
  ggplot(aes(x = year, y = width, col = type))+
  geom_line(linewidth = 1.1)+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
           length.out = 4) %>% 
         round()) 
  )+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("CI width"),
    x = "Year"
  )
ciWidthPltJrssb

fullCiPlotJrssb <- 
  ggarrange(
    ciPlotJrssb, ciWidthPltJrssb, common.legend = T, legend = "right", nrow = 2
  )
fullCiPlotJrssb
if(F){
  ggsave(plot = fullCiPlotJrssb, 
         filename = here("results/plotCaseStudyCbandsMainJrssb.pdf"),
         device = "pdf", width = 11, height = 5) #ursprl 4
}
ciTibPlt

#case study for the mean estimation of Precipitation----


##calculate estimators and bootstrap CIs
DbArrMean <- array( dim = c(nY - winSize + 1, 3))
SbArrMean <- array( dim = c(nY - winSize + 1, 3))
CbArrMean <- array( dim = c(nY - winSize + 1, 3))
bstVarsMean <- array(dim = c(nY - winSize + 1, 2))


if(F){#17mins with B = 10**3
  B <-  10**3
  t0 <- Sys.time()
  for(indW in seq(1, nY - winSize + 1)){
    dataWin <- (data %>% filter(
      year %in% seq(yearVec[indW], yearVec[indW] + 39)
    ))$RSK
    ciDbDatMean <- 
      ciCircmaxMean(dataWin, B = B, r = blockSize, k = 1, mthd = "db")
    DbArrMean[indW,c(1, 2,3)] <- ciDbDatMean #est, lower and upper
    
    
    ciCbDatMean <- 
      ciCircmaxMean(dataWin, B = B, r = blockSize, k = 2, mthd = "cb")
    CbArrMean[indW,c(1, 2,3)] <- ciCbDatMean #est, lower and upper
    
    if(indW %% 10 == 0) {print(indW)}
  }
  tDel <- difftime(Sys.time(), t0)
  cat(format(tDel))
  save(DbArrMean, CbArrMean, file = file.path("data", "array_mean_prcp_cs"))
}
load(file = file.path("data", "array_mean_cs_-0.10.4"))
# create tibble with data
ciTibMean <- bind_rows(
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "db",
    estim = DbArrMean[,1],
    lower = DbArrMean[,2],
    upper = DbArrMean[,3]
  ), 
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "cb",
    estim = CbArrMean[,1],
    lower = CbArrMean[,2],
    upper = CbArrMean[,3]
  ),
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "sb",
    estim = CbArrMean[,1],
    lower = SbArrMean[,2],
    upper = SbArrMean[,3]
  )
) %>% mutate(
  year = year + 40
)
ciTibMean %>% mutate(
  width = upper - lower
) 

avgWidthsMean <- (ciTibMean %>% mutate(
  width = upper - lower) %>% group_by(type) %>% 
    summarise(avgWidth = median(width)) %>% 
    filter(type != "sb") %>% select(avgWidth))[[1]]
avgWidthsMean[2]/avgWidthsMean[1]
##averaging of two year window
twoYeaAvgsMean <- 
  array(dim = c(length(ciTibMean$upper)/2, 2))
tmpTibMean <- tibble()
for(indT in seq(1,3)){
  indTT <- (ciTibMean$type %>% unique())[indT]
  updTibMean <- 
    (ciTibMean %>% filter(type == indTT))[,c(3,4,5)] %>% 
    apply(c(2), function(xx){stats::filter(xx, filter = c(1/2,1/2), sides = 1)}) %>% 
    as_tibble()
  updTibMean <- bind_cols(
    (ciTibMean %>% filter(type == indTT))[, c(1,2)],
    updTibMean
  )
  tmpTibMean <- tmpTibMean %>% bind_rows(updTibMean)
}
avgCiTibMean <- tmpTibMean %>% filter(year != 1919)
rm(tmpTibMean)



#plotting of confidence intervals

ciTibPltMean <- avgCiTibMean
##cb estim substituted by sb estimator
ciTibPltMean <- ciTibPltMean %>%  bind_cols(
  sbEstim = rep((ciTibPltMean %>% filter(type == "sb"))$estim, 3)
) %>% mutate(
  estim = ifelse(type == "cb", sbEstim, estim)
) %>% select(-c(sbEstim))

ciTibPltMean <- ciTibPltMean %>% filter(type != "sb")
ciTibPltMean$type <- 
  factor(ciTibPltMean$type, levels = c("db", "cb"))
ciTibPltMean$type <- 
  ciTibPltMean$type %>% factor(
    levels = c("cb", "db"),
    labels = c("cb" = "cb(2)", "db" = "db")
  )
ciTibPltMean

avgPrecDbMean <- (precDb %>% 
                stats::filter(sides = 1, filter = c(1/2,1/2)))[-1]

ciPlotMean <- 
  ciTibPltMean %>%  
  ggplot(aes(x = year))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.4,
              linewidth = 0.1, col = "black")+
  geom_line(
    data = ciTibPltMean,
    aes(x = year, y = estim, col = type), linewidth = 1)+
 # scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
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
ciPlotMean

ciTibPltMean <- ciTibPltMean %>% mutate(width = upper - lower) 
ciWidthPltMean <- 
  ciTibPltMean %>% 
  ggplot(aes(x = year, y = width, col = type))+
  geom_line(linewidth = 1.1)+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPltMean$year)[1], (ciTibPltMean$year)[nY - winSize], 
           length.out = 4) %>% 
         round()) 
  )+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("CI width"),
    x = "Year"
  )
ciWidthPltMean


fullCiPlotMean <- 
  ggarrange(
    ciPlotMean, ciWidthPltMean, common.legend = T, legend = "right"
  )
fullCiPlotMean

#explore end
if(F){
  ggsave(plot = fullCiPlotMean, filename = here("results/plotCaseStudyCbandsMain.pdf"),
         device = "pdf", width = 10, height = 5) #ursprl 4
}

#Jrssb

ciPlotMeanJrssb <- 
  ciTibPltMean %>% 
  ggplot(aes(x = year))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type, col = type), alpha = 0.4,
              linewidth = 0.5)+
  geom_line(
    data = ciTibPltMean,
    aes(x = year, y = estim, col = type), linewidth = 1.1)+
  #scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
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
ciPlotMeanJrssb

ciWidthPltMeanJrssb <- 
  ciTibPltMean %>% mutate(width = upper - lower) %>%  
  ggplot(aes(x = year, y = width, col = type))+
  geom_line(linewidth = 1.1)+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPltMean$year)[1], (ciTibPltMean$year)[nY - winSize], 
           length.out = 4) %>% 
         round()) 
  )+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("CI width"),
    x = "Year"
  )
ciWidthPltMeanJrssb

fullCiMeanPlotJrssb <- 
  ggarrange(
    ciPlotMeanJrssb, ciWidthPltMeanJrssb, common.legend = T, legend = "right", nrow = 2
  )
fullCiMeanPlotJrssb
if(F){
  ggsave(plot = fullCiMeanPlotJrssb, 
         filename = here("results/plotCaseStudyCbandsMeanMainJrssb.pdf"),
         device = "pdf", width = 11, height = 5) #ursprl 4
}

tib_mean_cs <- 
  bind_rows(ciTibPltMean %>% mutate(parameter = "mean"),
            ciTibPlt %>% mutate(parameter = "rl"))






library(ggplot2)
library(dplyr)

# 1️⃣ Daten für die zwei Facets anpassen
tibble_ci <- tib_mean_cs %>% mutate(panel = "Confidence Intervals")
tibble_ci$parameter <- factor(tibble_ci$parameter, levels = c("rl", "mean"), 
                              labels = c("Return Level", "Mean"))
tibble_width <- tib_mean_cs %>% mutate(panel = "CI Width")
tibble_width$parameter <- factor(tibble_width$parameter, levels = c("rl", "mean"), 
                                 labels = c("Return Level", "Mean"))
tibble_combined <- bind_rows(tibble_ci, tibble_width)


# 3️⃣ Plot mit `facet_wrap` und `scales = "free_y"`
plotCombCs <- 
  ggplot(tibble_combined) +
  
  # 🔹 1. Facet-Row: Estimate mit Konfidenzintervall
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = type), alpha = 0.4, data = tibble_ci, col = "black", linewidth = 0.1) +
  geom_line(aes(x = year, y = estim, color = type), linewidth = 1.1, 
            data = tibble_ci) +
  
  # 🔹 2. Facet-Row: CI Width
  geom_line(aes(x = year, y = width, color = type), linewidth = 1.1, 
            data = tibble_width) +
  
  # 🎨 Faceting nach `parameter` (Spalten) & `panel` (Zeilen)
  facet_wrap(vars(panel, parameter), scales = "free_y", ncol = length(unique(tibble_combined$parameter))) +
  
  # 🏷️ Labels & Design
  labs(
    x = "Year",
    y = "",
    color = "Bootstrap:",
    fill = "Bootstrap:",
    title = "Comparison of CI Estimation & Width"
  ) +
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
plotCombCs
if(F){
  ggsave(plot = plotCombCs, 
         filename = here("results/plotCaseStudy.pdf"),
         device = "pdf", width = 11, height = 7) #ursprl 4
}


#case study for the mean estimation of Temperature----


##calculate estimators and bootstrap CIs
DbArrMeanTmp <- array( dim = c(nY - winSize + 1, 3))
SbArrMeanTmp <- array( dim = c(nY - winSize + 1, 3))
CbArrMeanTmp <- array( dim = c(nY - winSize + 1, 3))
bstVarsMeanTmp <- array(dim = c(nY - winSize + 1, 2))

blockSizeTemp <- 92
if(F){#17mins with B = 10**3
  B <-  10**3
  t0 <- Sys.time()
  for(indW in seq(1, nY - winSize + 1)){
    dataWin <- (data %>% mutate(month = as.double(format(day, "%m"))) %>%
                  filter(year %in% seq(yearVec[indW], yearVec[indW] + 39), month %in% c(6,7,8)))$TXK
    ciDbDatMeanTmp <- 
      ciCircmaxMean(dataWin, B = B, r = blockSizeTemp, k = 1, mthd = "db")
    DbArrMeanTmp[indW,c(1, 2,3)] <- ciDbDatMeanTmp #est, lower and upper
    
    
    ciCbDatMeanTmp <- 
      ciCircmaxMean(dataWin, B = B, r = blockSizeTemp, k = 2, mthd = "cb")
    CbArrMeanTmp[indW,c(1, 2,3)] <- ciCbDatMeanTmp #est, lower and upper
    
    if(indW %% 10 == 0) {print(indW)}
  }
  tDel <- difftime(Sys.time(), t0)
  cat(format(tDel))
  save(DbArrMeanTmp, CbArrMeanTmp, file = file.path("data", "array_mean_tmp_cs"))
}
load(file = file.path("data", "array_mean_tmp_cs"))
# create tibble with data
ciTibMeanTmp <- bind_rows(
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "db",
    estim = DbArrMeanTmp[,1],
    lower = DbArrMeanTmp[,2],
    upper = DbArrMeanTmp[,3]
  ), 
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "cb",
    estim = CbArrMeanTmp[,1],
    lower = CbArrMeanTmp[,2],
    upper = CbArrMeanTmp[,3]
  ),
  tibble(
    year = yearVec[seq(1,nY - winSize + 1)],
    type = "sb",
    estim = CbArrMeanTmp[,1],
    lower = SbArrMeanTmp[,2],
    upper = SbArrMeanTmp[,3]
  )
) %>% mutate(
  year = year + 40
)
ciTibMeanTmp %>% mutate(
  width = upper - lower
) 

avgWidthsMeanTmp <- (ciTibMeanTmp %>% mutate(
  width = upper - lower) %>% group_by(type) %>% 
    summarise(avgWidth = median(width)) %>% 
    filter(type != "sb") %>% select(avgWidth))[[1]]
avgWidthsMeanTmp[2]/avgWidthsMeanTmp[1]
##averaging of two year window
twoYeaAvgsMeanTmp <- 
  array(dim = c(length(ciTibMeanTmp$upper)/2, 2))
tmpTibMeanTmp <- tibble()
for(indT in seq(1,3)){
  indTT <- (ciTibMeanTmp$type %>% unique())[indT]
  updTibMeanTmp <- 
    (ciTibMeanTmp %>% filter(type == indTT))[,c(3,4,5)] %>% 
    apply(c(2), function(xx){stats::filter(xx, filter = c(1/2,1/2), sides = 1)}) %>% 
    as_tibble()
  updTibMeanTmp <- bind_cols(
    (ciTibMeanTmp %>% filter(type == indTT))[, c(1,2)],
    updTibMeanTmp
  )
  tmpTibMeanTmp <- tmpTibMeanTmp %>% bind_rows(updTibMeanTmp)
}
avgCiTibMeanTmp <- tmpTibMeanTmp %>% filter(year != 1919)
rm(tmpTibMeanTmp)



#plotting of confidence intervals

ciTibPltMeanTmp <- avgCiTibMeanTmp
##cb estim substituted by sb estimator
ciTibPltMeanTmp <- ciTibPltMeanTmp %>%  bind_cols(
  sbEstim = rep((ciTibPltMeanTmp %>% filter(type == "sb"))$estim, 3)
) %>% mutate(
  estim = ifelse(type == "cb", sbEstim, estim)
) %>% select(-c(sbEstim))

ciTibPltMeanTmp <- ciTibPltMeanTmp %>% filter(type != "sb")
ciTibPltMeanTmp$type <- 
  factor(ciTibPltMeanTmp$type, levels = c("db", "cb"))
ciTibPltMeanTmp$type <- 
  ciTibPltMeanTmp$type %>% factor(
    levels = c("cb", "db"),
    labels = c("cb" = "cb(2)", "db" = "db")
  )
ciTibPltMeanTmp

avgPrecDbMeanTmp <- (precDb %>% 
                    stats::filter(sides = 1, filter = c(1/2,1/2)))[-1]

ciPlotMeanTmp <- 
  ciTibPltMeanTmp %>%  
  ggplot(aes(x = year))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.4,
              linewidth = 0.1, col = "black")+
  geom_line(
    data = ciTibPltMeanTmp,
    aes(x = year, y = estim, col = type), linewidth = 1)+
  # scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
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
ciPlotMeanTmp

ciTibPltMeanTmp <- ciTibPltMeanTmp %>% mutate(width = upper - lower) 
ciWidthPltMeanTmp <- 
  ciTibPltMeanTmp %>% 
  ggplot(aes(x = year, y = width, col = type))+
  geom_line(linewidth = 1.1)+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPltMeanTmp$year)[1], (ciTibPltMeanTmp$year)[nY - winSize], 
           length.out = 4) %>% 
         round()) 
  )+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("CI width"),
    x = "Year"
  )
ciWidthPltMeanTmp


fullCiPlotMeanTmp <- 
  ggarrange(
    ciPlotMeanTmp, ciWidthPltMeanTmp, common.legend = T, legend = "right"
  )
fullCiPlotMeanTmp

#explore end
if(F){
  ggsave(plot = fullCiPlotMeanTmp, filename = here("results/plotCaseStudyCbandsMain.pdf"),
         device = "pdf", width = 10, height = 5) #ursprl 4
}

#Jrssb

ciPlotMeanTmpJrssb <- 
  ciTibPltMeanTmp %>% 
  ggplot(aes(x = year))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type, col = type), alpha = 0.4,
              linewidth = 0.5)+
  geom_line(
    data = ciTibPltMeanTmp,
    aes(x = year, y = estim, col = type), linewidth = 1.1)+
  #scale_y_continuous(limits = c(0,70))+
  scale_x_continuous(
    breaks = 
      seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
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
ciPlotMeanTmpJrssb

ciWidthPltMeanTmpJrssb <- 
  ciTibPltMeanTmp %>% mutate(width = upper - lower) %>%  
  ggplot(aes(x = year, y = width, col = type))+
  geom_line(linewidth = 1.1)+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPltMeanTmp$year)[1], (ciTibPltMeanTmp$year)[nY - winSize], 
           length.out = 4) %>% 
         round()) 
  )+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("CI width"),
    x = "Year"
  )
ciWidthPltMeanTmpJrssb

fullCiMeanTmpPlotJrssb <- 
  ggarrange(
    ciPlotMeanTmpJrssb, ciWidthPltMeanTmpJrssb, common.legend = T, legend = "right", nrow = 2
  )
fullCiMeanTmpPlotJrssb
if(F){
  ggsave(plot = fullCiMeanTmpPlotJrssb, 
         filename = here("results/plotCaseStudyCbandsMeanTmpMainJrssb.pdf"),
         device = "pdf", width = 11, height = 5) #ursprl 4
}

tib_mean_cs <- 
  bind_rows(ciTibPltMeanTmp %>% mutate(parameter = "mean"),
            ciTibPlt %>% mutate(parameter = "rl"))






library(ggplot2)
library(dplyr)

# 1️⃣ Daten für die zwei Facets anpassen
tibble_ci <- tib_mean_cs %>% mutate(panel = "Confidence Intervals")
tibble_ci$parameter <- factor(tibble_ci$parameter, levels = c("rl", "mean"), 
                              labels = c("Return Level", "MeanTmp"))
tibble_width <- tib_mean_cs %>% mutate(panel = "CI Width")
tibble_width$parameter <- factor(tibble_width$parameter, levels = c("rl", "mean"), 
                                 labels = c("Return Level", "MeanTmp"))
tibble_combined <- bind_rows(tibble_ci, tibble_width)


# 3️⃣ Plot mit `facet_wrap` und `scales = "free_y"`
plotCombCs <- 
  ggplot(tibble_combined) +
  
  # 🔹 1. Facet-Row: Estimate mit Konfidenzintervall
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = type), alpha = 0.4, data = tibble_ci, col = "black", linewidth = 0.1) +
  geom_line(aes(x = year, y = estim, color = type), linewidth = 1.1, 
            data = tibble_ci) +
  
  # 🔹 2. Facet-Row: CI Width
  geom_line(aes(x = year, y = width, color = type), linewidth = 1.1, 
            data = tibble_width) +
  
  # 🎨 Faceting nach `parameter` (Spalten) & `panel` (Zeilen)
  facet_wrap(vars(panel, parameter), scales = "free_y", ncol = length(unique(tibble_combined$parameter))) +
  
  # 🏷️ Labels & Design
  labs(
    x = "Year",
    y = "",
    color = "Bootstrap:",
    fill = "Bootstrap:",
    title = "Comparison of CI Estimation & Width"
  ) +
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
plotCombCs
if(F){
  ggsave(plot = plotCombCs, 
         filename = here("results/plotCaseStudy.pdf"),
         device = "pdf", width = 11, height = 8) #ursprl 4
}