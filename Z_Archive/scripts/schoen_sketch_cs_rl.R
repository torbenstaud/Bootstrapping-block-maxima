#load libraries----
library(readr) #read and download files
library(tidyr) #data wrangling
library(dplyr) #data wrangling
library(ggplot2) #plotting
library(ggpubr)#arranging 
library(here) #relative paths
library(Rcpp) #cpp api
library(magrittr)
#source files----
source(here("src/0source.R"))
source(here("src/0source_ts.R"))
sourceCpp(here("src/0sourceC.cpp"))
sourceCpp(here("src/0sourceCn.cpp"))
source(here("src/sourceCase.R"))
source(here("src/theme.R"))




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
  c("sb-cb" = "#F8766D",
    "cb(2)" = "#F8766D",  
    "cb(3)" = "#7CAE00",  
    "db" = "#00BFC4",  
    "sb" = "#C77CFF")





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
  level <- 0.05
  for(indW in seq(1, nY - winSize + 1)){
    dataWin <- (data %>% filter(
      year %in% seq(yearVec[indW], yearVec[indW] + 39)
    ))$RSK
    ciDbDat <- 
      ciCircmaxRl(dataWin, B = B, r = blockSize, k = 1, mthd = "db", onesided = F,
                  botcutoffGam = -0.3, topcutoffGam = 0.5, niv = level)
    # ciCircmaxRlPos(dataWin, B = B, r = blockSize, k = 1, mthd = "db", onesided = F,
    #             botcutoffGam = -0.3, topcutoffGam = 0.5, niv = level) #werte für die cutoffs, da die shape über den ganzen Zeitraum als ~0.1 geschätzt wird und wir ein 0.4 band drum legen
    DbArr[indW,c(1, 2,3)] <- ciDbDat #est, lower and upper
    
    if(max(abs(ciDbDat)) >= 100 ){
      sprintf("Error at indW = %g für db", indW)
    }
    
    #cat("db fertig") #debug
    ciCbDat <- 
      ciCircmaxRl(dataWin, B = B, r = blockSize, k = 2, mthd = "cb", onesided = F, 
                  botcutoffGam = -0.3, topcutoffGam = 0.5,
                  niv = level)
    # ciCircmaxRlPos(dataWin, B = B, r = blockSize, k = 2, mthd = "cb", onesided = F, 
    #                botcutoffGam = -0.3, topcutoffGam = 0.5,
    #             niv = level)
    CbArr[indW,c(1, 2,3)] <- ciCbDat #est, lower and upper
    
    if(max(abs(ciCbDat)) >= 100 ){
      sprintf("Error at indW = %g für cb", indW)
    }
    
    if(indW %% 10 == 0) {print(indW)}
  }
  tDel <- difftime(Sys.time(), t0)
  cat(format(tDel))
  #save(DbArr, CbArr, file = file.path("data", "array_rl_cs-0.30.5Niv0.9"))
}
load(file = file.path("data", "array_rl_cs-0.30.5")) # DbArr, CbArr
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
ciTib
##corecct the boundaries
load("data/coeff_rl0.5")
shape_est <- (precDb %>% fgev())$par[3]
corr_factor <- coeff_rl[1]+shape_est*coeff_rl[2]+ winSize *coeff_rl[3]

ciTib <- 
  ciTib %>% mutate(
    lower = (lower - estim)*corr_factor + lower, 
    upper = (upper - estim)*corr_factor + upper
  ) 

avgWidths <- (ciTib %>% mutate(
  width = upper - lower) %>% group_by(type) %>% 
    summarise(avgWidth = median(width)) %>% 
    filter(type != "sb") %>% select(avgWidth))[[1]]
avgWidths[1]/avgWidths[2] #wenn < 1 gut!
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



yearBounds <- 
  (ciTibPlt$year)[c(1, nY - winSize)]



#case study for the mean estimation of Precipitation----


##calculate estimators and bootstrap CIs
DbArrMean <- array( dim = c(nY - winSize + 1, 3))
SbArrMean <- array( dim = c(nY - winSize + 1, 3))
CbArrMean <- array( dim = c(nY - winSize + 1, 3))
bstVarsMean <- array(dim = c(nY - winSize + 1, 2))


if(F){#17mins with B = 10**3
  B <-  10**2
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
load(file = file.path("data", "array_mean_prcp_cs"))
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
avgWidthsMean[1]/avgWidthsMean[2] #wenn < 1: cb besser als db
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

##corecct the boundaries
load("data/coeff_rl")
corr_factor_mean <- coeff_mean[1]+shape_est*coeff_mean[2]+ winSize *coeff_mean[3]

#plotting of confidence intervals

ciTibPltMean <- avgCiTibMean %>% 
  mutate(lower = (lower -estim)*corr_factor_mean + lower, 
         upper = (upper -estim)*corr_factor_mean+ upper)
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


tib_mean_cs <- 
  bind_rows(ciTibPltMean %>% mutate(parameter = "mean"),
            ciTibPlt %>% mutate(parameter = "rl"))








# 1️⃣ Daten für die zwei Facets anpassen
tibble_ci <- tib_mean_cs %>% mutate(panel = "Confidence Interval")
tibble_ci$parameter <- factor(tibble_ci$parameter, levels = c("rl", "mean"), 
                              labels = c("Return Level", "Mean"))
tibble_width <- tib_mean_cs %>% mutate(panel = "CI Width")
tibble_width$parameter <- factor(tibble_width$parameter, levels = c("rl", "mean"), 
                                 labels = c("Return Level", "Mean"))
type_levels <- c("db", "cb(2)")
type_labels <- c("db" = "db", "cb(2)" = "sb-cb")

# In ALLEN drei Dataframes:
set_panel_order <- function(df) {
  df$panel <- factor(as.character(df$panel),
                     levels = c("Confidence Interval", "CI Width"))
  df$type <- factor(as.character(df$type),
                    levels = type_levels,
                    labels = type_labels)
  df
}

tibble_ci       <- set_panel_order(tibble_ci)
tibble_width    <- set_panel_order(tibble_width) %>% 
  mutate(width = upper - lower)

tibble_combined <- 
  bind_rows(tibble_ci, tibble_width)




plotRlCsPrcp <- 
  ggplot(tibble_combined %>% filter(parameter == "Return Level")) +
  
  # 🔹 1. Facet-Row: Estimate mit Konfidenzintervall
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = type), alpha = 0.4, data = tibble_ci %>% filter(parameter == "Return Level"), col = "black", linewidth = 0.1) +
  geom_line(aes(x = year, y = estim, color = type), linewidth = 1.1, 
            data = tibble_ci %>% filter(parameter == "Return Level")) +
  
  # 🔹 2. Facet-Row: CI Width
  geom_line(aes(x = year, y = width, color = type), linewidth = 1.1, 
            data = tibble_width %>% filter(parameter == "Return Level")) +
  
  # 🎨 Faceting nach `parameter` (Spalten) & `panel` (Zeilen)
  facet_grid(rows = vars(panel), cols = vars(parameter), scales = "free_y") +
  
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background.y = element_blank(),
    strip.text.y = element_blank()
  )+
  scale_color_manual(
    values = ownPalette,
    labels = type_labels
  ) +
  scale_fill_manual(
    values = ownPalette,
    labels = type_labels
  )


plotRlCsPrcp
#zweiter Plot 2x1 Mean
plotMeanCsPrcp <- 
  ggplot(tibble_combined %>% filter(parameter == "Mean")) +
  
  # 🔹 1. Facet-Row: Estimate mit Konfidenzintervall
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = type), alpha = 0.4, data = tibble_ci %>% filter(parameter == "Mean"), col = "black", linewidth = 0.1) +
  geom_line(aes(x = year, y = estim, color = type), linewidth = 1.1, 
            data = tibble_ci %>% filter(parameter == "Mean")) +
  
  # 🔹 2. Facet-Row: CI Width
  geom_line(aes(x = year, y = width, color = type), linewidth = 1.1, 
            data = tibble_width %>% filter(parameter == "Mean")) +
  
  # 🎨 Faceting nach `parameter` (Spalten) & `panel` (Zeilen)
  facet_grid(rows = vars(panel), cols = vars(parameter), scales = "free_y") +
  
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
  )+
  scale_color_manual(
    values = ownPalette,
    labels = type_labels
  ) +
  scale_fill_manual(
    values = ownPalette,
    labels = type_labels
  )

plotMeanCsPrcp

plot_combinedCsFaCoRlMean <- plotRlCsPrcp + plotMeanCsPrcp +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "right")

plot_combinedCsFaCoRlMean


if(F){
  ggsave(plot = plot_combinedCsFaCoRlMean, 
         filename = here("results/plot_combinedCsFaCoRlMean.pdf"),
         device = "pdf", width = 13, height = 7) #ursprl 4
}
