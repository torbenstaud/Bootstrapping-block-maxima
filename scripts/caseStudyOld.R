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
data <- dataRaw %>% select(MESS_DATUM, RSK) %>% 
  rename(day = "MESS_DATUM") %>% 
  mutate(day = as.Date(as.character(day), format = "%Y%m%d"))
#rm(dataRaw)
data
## missing or absurd data?
data %>% na.omit()
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
data <- data %>% filter(RSK >= 0) #starts at 1879
data %>% ggplot(aes(x = day, y = RSK))+
  geom_line()
data <- data %>%  mutate(year = as.double(format(day, "%Y")))
nY <- max(data$year) - min(data$year) + 1
## -> left with 145 full years of daily observations (2023-1879+1)

## diagnostics: are disjoint block maxima over a year approx independent?
dbmsFull <- data %>% group_by(year) %>% summarize(blockmax = max(RSK))
dbmsFull$blockmax %>% acf() #no significant evidence against independence
#rolling window over 40 years and 1 year jumps
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
) %>% mutate(
  year = year + 40
)

avgWidths <- (ciTib %>% mutate(
  width = upper - lower) %>% group_by(type) %>% 
  summarise(avgWidth = median(width)) %>% 
  filter(type != "sb") %>% select(avgWidth))[[1]]
#avgWidths[2]/avgWidths[1] 7% smaller ci on average
##filtering (size 2, weights unif)
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



##filtered (window 2, weights: unif) db block maxima
avgPrecDb <- (precDb %>% 
  stats::filter(sides = 1, filter = c(1/2,1/2)))[-1]

##plotting
ciPlot <- 
  ciTibPlt %>% 
  ggplot(aes(x = year))+
  #geom_line(
  #  data = tibble(
  #    x = seq(yearBounds[1], yearBounds[2]),
  #    y = avgPrecDb[seq(1, yearBounds[2] - yearBounds[1] +1)]
  #  ),
  #  aes(x = x, y = (y - min(avgPrecDb))/5), col = "black"
  #)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type, col = type), alpha = 0.4,
              linewidth = 0.1)+
  geom_line(
    data = ciTibPlt,
    aes(x = year, y = estim, col = type), linewidth = 1)+
  scale_y_continuous(limits = c(0,70))+
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


## plot filtered(moving average of window 2 with weights 1/2) CI widths
ciTibPlt <- ciTibPlt %>% mutate(width = upper - lower) 
ciWidthPlt <- 
ciTibPlt %>% 
  ggplot(aes(x = year, y = width, col = type))+
  geom_line()+
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

##full plot for CIs and CI width
fullCiPlot <- 
  ggarrange(
    ciPlot, ciWidthPlt, common.legend = T, legend = "right"
  )
fullCiPlot
if(F){
  ggsave(plot = fullCiPlot, filename = here("results/plotCaseStudyCbandsMain.pdf"),
         device = "pdf", width = 10, height = 5) #ursprl 4
}
##relative CI width
relCiTib <- ciTibPlt %>% filter(type == "cb(2)") %>% bind_cols(
  tibble( dbWidth = (ciTibPlt %>% filter(type == "db"))$width)
) %>% mutate(relWidth = dbWidth/width) %>% 
  select(c(year, relWidth))
relCiTib

ciRelWidthPlt <- 
  relCiTib %>% 
  ggplot(aes(x = year, y = relWidth))+
  geom_line(col = "#F8766D")+
  scale_x_continuous(
    breaks = 
      (seq((ciTibPlt$year)[1], (ciTibPlt$year)[nY - winSize], 
           length.out = 4) %>% 
         round()) 
  )+
  scale_y_continuous(
    limits = c(1, 1.4)
  )+
  geom_hline(
    yintercept = 1, linetype = "longdash", col = "#00BFC4"
  )+
  themePlot+
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  labs(
    title = "Precipitation heigth in ",
    y = paste0("Relative CI width"),
    x = "Year"
  )
ciRelWidthPlt

##full plot for CIs and RELATIVE CI width (both filtered)
fullCiPlot2 <- 
  ggarrange(
    ciPlot, ciRelWidthPlt, common.legend = T, legend = "right"
  )
fullCiPlot2

if(F){
  ggsave(plot = fullCiPlot2, filename = here("results/plotCaseStudyCbandsRel.pdf"),
         device = "pdf", width = 10, height = 5) #ursprl 4
}



