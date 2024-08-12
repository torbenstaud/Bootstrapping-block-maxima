# --------------------------------------------------------------------
# Script Name: speedComparison.R
# Purpose: This scripts produces runtime comparison plots for disjoint blocks 
#          vs circmax bootstrapping.
# Author: Torben Staud
# Date: 2024-08
# --------------------------------------------------------------------

# Clear the environment and console
rm(list=ls())  # Clear all objects from the current R environment
cat("\014")    # Clear the console

# ==== SECTION 1: Load required libraries ======================================
# Description: Load all the R packages required for the script.
#_______________________________________________________________________________
library(quantmod) #download financial data from yahoo finance
library(tidyverse) #tidyverse package suite for data manipulation
require(here) #construct paths that are relative to the project's root directory
library(Rcpp) #enable the use of C++ code within R for performance
library(gsl) # package for zeta functions etc.


# ==== SECTION 2: SOURCE EXTERNAL FILES ===========================================
# Description: Source external scripts which might include custom functions
#              or specific project configurations.
#______________________________________________________________________________
source(here("src/0source.r"))
sourceCpp(here("src/0sourceC.cpp"))
source(here("src/sourceCase.R"))
source(here("src/theme.R"))

# ==== SECTION 3: main script ===========================================
# Description: code for the case study data and plots
#______________________________________________________________________________
start_date <- as.Date("1974-01-01")
end_date <- as.Date("2023-12-31")
data0 <- getSymbols("^GSPC", src = "yahoo", from = start_date, 
                    to = end_date, auto.assign = FALSE)
data0_tibble <- as_tibble(data.frame(Date = index(data0), coredata(data0)))
logrets <- data0_tibble %>% mutate(log_return = log(GSPC.Close/lag(GSPC.Close))) %>% 
  drop_na(log_return) %>% mutate(log_wins = pmax(log_return,0), log_loss = -pmin(log_return,0)) %>% 
  select(c(Date, log_wins, log_loss)) 
plogrets <- logrets$log_wins
nlogrets <- logrets$log_loss
# Sample Size:
n <- length(logrets$log_wins)
# Dates:
date_all <- as.Date(index(data0))
date_years <- format(date_all, "%Y")
date_quarters <- quarters(date_all)
date_months <- months(date_all)
date_final <- paste(date_years,":",date_quarters, sep="")
date_final_m <- paste(date_years,":",date_months, sep="")  


# Number of quarters:
l.q <- 50*4 #50 years correspond to 50*4 quarters

# Calculate Index where quarter changes (variable qc), and the quarter number (from 1 to l.q=200) 
# each observation belongs to (variable date.quarters2)
qc <- numeric(l.q)
date_quarters2 <- numeric(n)
date_quarters2[1] <- 1; qc[1] <- 1  # Initialize
z <- 1
for(i in seq(2,n)) {
  if(date_quarters[i] != date_quarters[i-1]) {
    z <- z+1
    qc[z] <- i
  }
  date_quarters2[i] <- z
}
# Calculate Quarterly Maxima (disjoint):
plogrets_dbm <- numeric(l.q)
nlogrets_dbm <- numeric(l.q)
for(i in 1:(l.q-1)) {
  plogrets_dbm[i] <- max(plogrets[qc[i]:(qc[i+1]-1)])
  nlogrets_dbm[i] <- max(nlogrets[qc[i]:(qc[i+1]-1)])
}
plogrets_dbm[l.q] <- max(plogrets[qc[l.q]])
nlogrets_dbm[l.q] <- max(nlogrets[qc[l.q]])


#initialize estimation arrays
pml.db <- matrix(nrow=l.q - 40, ncol=2)
nml.db <- matrix(nrow=l.q - 40, ncol=2)
pml.sb <- matrix(nrow=l.q - 40, ncol=2)
nml.sb <- matrix(nrow=l.q - 40, ncol=2)
pCiUpper <- matrix(nrow=l.q - 40, ncol=2)
pCiLower <- matrix(nrow=l.q - 40, ncol=2)
nCiUpper <- matrix(nrow=l.q - 40, ncol=2)
nCiLower <- matrix(nrow=l.q - 40, ncol=2)
pVars <- matrix(nrow = l.q - 40, ncol = 2)
nVars <- matrix(nrow = l.q - 40, ncol = 2)

# length of quarter: (62 rougly corresonds to one trading quarter)
r <- 62 
# Estimation loop (rolling window over the quarters)
if(F){
  
for(i in 1:(l.q - 40)) {
  ind <- which( date_quarters2 %in% (i:(i+39)) )
  pdata <- plogrets[ind]; 
  ndata <- nlogrets[ind]
  # Disjoint Blocks:
  ##To do für mich:
  # -Eigene Mle und slidmax funktionen drauf anwenden
  # -bei slidmax loopen verhindern, -bei mle drauf achten, dass alles passt
  ##
  pml.db[i,] <- mleFre(plogrets_dbm[i:(i+39)])
  nml.db[i,] <- mleFre(nlogrets_dbm[i:(i+39)])  
  # Sliding Blocks:
  pml.sb[i,] <- mleFre(slidMaxCNoLoop(pdata,r))
  nml.sb[i,] <- mleFre(slidMaxCNoLoop(ndata,r))
  pQuantsAndVars <- pdata %>% slidMaxCNoLoop(r) %>% ciCircmax(r,2, B = 2*10^3)
  nQuantsAndVars <- ndata %>% slidMaxCNoLoop(r) %>% ciCircmax(r,2, B = 2*10^3)
  pCiBst <- pQuantsAndVars[[1]]
  nCiBst <- nQuantsAndVars[[1]] 
  pVars[i,] <- pQuantsAndVars[[2]]
  nVars[i,] <- nQuantsAndVars[[2]]
  pCiUpper[i,] <- pCiBst[,2] #scale cis win
  pCiLower[i,] <- pCiBst[,1] #shape cis win
  nCiUpper[i,] <- nCiBst[,2] #scale cis loss
  nCiLower[i,] <- nCiBst[,1] #shape cis loss
}
  save(
    pml.db, nml.db, pml.sb, nml.sb, pCiUpper, pCiLower, nCiUpper, nCiLower, 
    pVars, nVars, 
    file = here("data/caseStudyData.R")
  )
}
load(here("data/caseStudyData.R"))

nDbTib <- tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  method = "db",
  type = "losses",
  shape = nml.db[,1],
  scale = nml.db[,2]
)
pDbTib <- tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  method = "db",
  type = "wins",
  shape = pml.db[,1],
  scale = pml.db[,2]
)
nSlTib <- tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  method = "sb",
  type = "losses",
  shape = nml.sb[,1],
  scale = nml.sb[,2]
)
pSlTib <- tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  method = "sb",
  type = "wins",
  shape = pml.sb[,1],
  scale = pml.sb[,2]
)
#we discard db from now on
#create bootstrap confidence interval tibbles
pCiTib <- bind_rows(tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  type = "wins",
  param = "scale",
  ciType = "bstr",
  Upper = pCiUpper[,2],
  Lower = pCiUpper[,1]
), tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  type = "wins",
  param = "shape",
  ciType = "bstr",
  Upper = pCiLower[,2],
  Lower = pCiLower[,1]
))  %>% left_join(
  pSlTib %>% pivot_longer(cols = c(shape, scale), names_to = "param",
                          values_to = "val")) %>% 
  select(-method)

nCiTib <- bind_rows(tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  type = "losses",
  param = "scale",
  ciType = "bstr",
  Upper = nCiUpper[,2],
  Lower = nCiUpper[,1]
), tibble(
  data = nlogrets_dbm[41:200],
  dateInd = seq_along(date_final[qc[41:200]]),
  date = date_final[qc[41:200]],
  type = "losses",
  param = "shape",
  ciType = "bstr",
  Upper = nCiLower[,2],
  Lower = nCiLower[,1]
))  %>% left_join(
  nSlTib %>% pivot_longer(cols = c(shape, scale), names_to = "param",
                          values_to = "val")) %>% 
  select(-method)
##CIs based on normal approximation with true asymptotic variance formula
pCiAsyTib <- 
  pSlTib %>% 
  mutate(asyQuant_Shape = qnorm(0.975, mean=0, 
                                sd = sqrt(v11(shape)/ 40)),
         asyQuant_Scale = qnorm(0.975, mean=0, 
                                sd = scale * sqrt(v22(shape)/ 40))
  ) %>% 
  pivot_longer(cols = c(shape, scale), values_to = "val", names_to = "param") %>% 
  rename(shape = "asyQuant_Shape", scale = "asyQuant_Scale") %>% 
  pivot_longer(cols = c(shape, scale), values_to = "asyQuant") %>%
  mutate(ciType = "nApprox", Upper = val + asyQuant, Lower = val - asyQuant) %>% 
  filter(param == name) %>% select(-c(name, method, asyQuant))
nCiAsyTib <- 
  nSlTib %>% 
  mutate(asyQuant_Shape = qnorm(0.975, mean=0, 
                                sd = sqrt(v11(shape)/ 40)),
         asyQuant_Scale = qnorm(0.975, mean=0, 
                                sd = scale * sqrt(v22(shape)/ 40))
  ) %>% 
  pivot_longer(cols = c(shape, scale), values_to = "val", names_to = "param") %>% 
  rename(shape = "asyQuant_Shape", scale = "asyQuant_Scale") %>% 
  pivot_longer(cols = c(shape, scale), values_to = "asyQuant") %>%
  mutate(ciType = "nApprox", Upper = val + asyQuant, Lower = val - asyQuant) %>% 
  filter(param == name) %>% select(-c(name, method, asyQuant))
##CIs based on normal approx with variance bootstrapped

pCiAsyBstrTib <- 
  pSlTib %>% bind_cols(tibble(shapeV = pVars[,1], scaleV = pVars[,2])) %>% 
  mutate(asyQuant_Shape = qnorm(0.975, mean=0, 
                                sd = sqrt(shapeV)),
         asyQuant_Scale = qnorm(0.975, mean=0, 
                                sd = sqrt(scaleV))
  ) %>% 
  pivot_longer(cols = c(shape, scale), values_to = "val", names_to = "param") %>% 
  rename(shape = "asyQuant_Shape", scale = "asyQuant_Scale") %>% 
  pivot_longer(cols = c(shape, scale), values_to = "asyQuant") %>%
  mutate(ciType = "bstrApprox", Upper = val + asyQuant, Lower = val - asyQuant) %>% 
  filter(param == name) %>% select(-c(name, method, asyQuant))
nCiAsyBstrTib <- 
  nSlTib %>% bind_cols(tibble(shapeV = nVars[,1], scaleV = nVars[,2])) %>% 
  mutate(asyQuant_Shape = qnorm(0.975, mean=0, 
                                sd = sqrt(shapeV)),
         asyQuant_Scale = qnorm(0.975, mean=0, 
                                sd = sqrt(scaleV))
  ) %>% 
  pivot_longer(cols = c(shape, scale), values_to = "val", names_to = "param") %>% 
  rename(shape = "asyQuant_Shape", scale = "asyQuant_Scale") %>% 
  pivot_longer(cols = c(shape, scale), values_to = "asyQuant") %>%
  mutate(ciType = "bstrApprox", Upper = val + asyQuant, Lower = val - asyQuant) %>% 
  filter(param == name) %>% select(-c(name, method, asyQuant))


fullCiTib <- bind_rows(
  bind_rows(pCiTib, nCiTib), bind_rows(pCiAsyTib, nCiAsyTib),
  bind_rows(pCiAsyBstrTib, nCiAsyBstrTib)
)


rescale <- 100

plotTib <- fullCiTib %>% 
  mutate(val = ifelse(param == "scale", val*rescale, val),
         Upper = ifelse(param == "scale", Upper*rescale, Upper), 
         Lower = ifelse(param == "scale", Lower*rescale, Lower))


###Plotting
#textsize = 15 from theme.R
cBandsPlot <- 
plotTib %>% 
  filter(param == "shape") %>% 
  ggplot(aes(x = dateInd, val, col = ciType))+
  geom_line(col = "black")+
  geom_line(data= plotTib, aes(x = dateInd, y = data*5), col = "darkblue")+
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2,
              linewidth = 0.9)+
  facet_wrap(~type)+
  scale_x_continuous(breaks = c(1,49,97,145), 
                     labels = (plotTib$date %>% unique)[c(1,49,97,145)])+
  facet_wrap(~type, 
             labeller = labeller(type = c(losses = "Losses", wins = "Wins")))+
  scale_color_discrete(name = "Type", 
                       labels = 
                         c(bstr = "BCI", bstrApprox = "NBCI", nApprox = "NCI")
  )+
  themePlot+
  theme(
    axis.title.x = element_blank()
  )+
  labs(
    title = "S&P 500 shape parameter confidence bands",
    y = "Shape"
  )
cBandsPlot

if(F){
ggsave(plot = cBandsPlot, filename = here("results/plotCaseStudyCbands.pdf"),
       device = "pdf", width = 12, height = 5)
}


##Plot of the widths
plotTib %>% 
  filter(param == "shape") %>% 
  mutate(width = Upper - Lower) %>% 
  ggplot(aes(x = dateInd, width, col = ciType))+
  geom_line()+
  facet_wrap(~type, 
             labeller = labeller(type = c(losses = "Losses", wins = "Wins")))+
  scale_x_continuous(breaks = c(1,49,97,145), 
                     labels = (plotTib$date %>% unique)[c(1,49,97,145)])+
  scale_color_discrete(name = "Type", 
                       labels = 
                         c(bstr = "BCI", bstrApprox = "NBCI", nApprox = "NCI")
  )

##average widths
plotTib %>% group_by(ciType, type) %>% 
  filter(param == "shape") %>% 
  mutate(width = Upper - Lower) %>% 
  summarise(avgWidth = mean(width)) %>% 
  ggplot(aes(x = ciType, y = avgWidth))+
  geom_col()+
  facet_wrap(~type, 
             labeller = labeller(type = c(losses = "Losses", wins = "Wins")))+
  themePlot+
  labs(y = "Average width of confidence interval")+
  scale_x_discrete(name = "Method",
                   labels = 
                     c(bstr = "BCI", bstrApprox = "NBCI", nApprox = "NCI")
    )

