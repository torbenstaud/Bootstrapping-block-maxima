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
require(Rcpp) #enable the use of C++ code within R for performance
require(here) #construct paths that are relative to the project's root directory
require(microbenchmark) #accurate timing of code execution
require(ggthemes) #for a colorblind-friendly palette
require(tidyverse) #tidyverse package suite for data manipulation
require(ggpubr) #arranging multiple plots

# ==== SECTION 2: SOURCE EXTERNAL FILES ===========================================
# Description: Source external scripts which might include custom functions
#              or specific project configurations.
#______________________________________________________________________________

source(here("src/theme.R")) #plot theme


# ==== SECTION 3: MAIN ANALYSIS ================================================
# Description: The core analytical processes, including data loading, cleaning,
#              processing, analysis, and visualization.
#______________________________________________________________________________


## ==== SECTION 3.1: fixed r (block size) ======================================
# Description: The core analytical processes, including data loading, cleaning,
#              processing, analysis, and visualization.
#_______________________________________________________________________________

#load the processed data file; name: fullBenchTibR
load(here("data/benchTibMleRFix"))


labsPlot <- labs(
  x = "Effective sample size m",
  col = "Bootstrap")
###-factors  ========================================================
tmpDbTibMle <- fullBenchTibR %>% filter(expr == "disjoint") %>% rename(timeDB = "timeMed") %>% 
  select(-expr)
benchTibRelMle <- left_join(fullBenchTibR, tmpDbTibMle) %>% filter(expr !="disjoint") %>% 
  mutate(relTime = timeMed/timeDB)
benchTibRelMle$expr <- benchTibRelMle$expr %>% factor(
  levels = c("sliding", "k2", "k3",  "k2No", "k3No","disjoint"), 
  labels = c("sb", "cb(2)", "cb(3)", "k2No", "k3No", "db")
)
ownPalette <- #based on dark2
  c("cb(2)" = "#F8766D",  
    "cb(3)" = "#7CAE00",  
    "db" = "#00BFC4",  
    "sb" = "#C77CFF")

benchTibRelMle$B <- benchTibRelMle$B %>% factor(
  levels = c(250, 500, 750, 1000), 
  labels = c("B = 250", "B = 500", "B = 750", "B = 1000")
)
textSize <- 17
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
###plotting----
relTimePlot <-
benchTibRelMle  %>% mutate(n = m*r) %>% 
  filter(!expr %in% c("k2No", "k3No"), m >= 40) %>% 
  ggplot(aes(x = n, y = relTime, col = expr))+ 
  geom_line(linewidth = 1.1)+
  geom_hline(yintercept = 1, col = "#00BFC4", linetype = "longdash")+
  facet_wrap(~B, ncol = 4)+
  labsPlot+
  labs(x = "Total sample size n",
       title = "Relative runtimes for Fréchet MLE", 
       y = "Relative time")+
  themePlot+
  scale_color_manual(values = ownPalette)+
  scale_x_continuous(
    breaks = c(4000, 6000, 8000)
  )+
  theme( 
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  )
relTimePlot
if(F){
ggsave(relTimePlot, path = here("results/"), 
       filename = "plotrelTimePlotMleRFix.pdf", device = "pdf",
       width = 10, height = 5, 
)
}

## ==== SECTION 3.1: fixed n (sample length) ======================================
# Description: The core analytical processes, including data loading, cleaning,
#              processing, analysis, and visualization.
#_______________________________________________________________________________

load(here("data/benchTibMleNFix")) #name: fullBenchTibN
fullBenchTibN %>% 
  filter(!expr %in% c("k2No", "k3No")) %>% 
  ggplot(aes(x = m, y = timeMed/1000, col = expr))+ geom_line()+
  facet_wrap(~B, scales = "free_y")+ 
  labs(x = "Effective sample size (number of disjoint blocks)", 
       y = "Time in seconds",
       title = "Median Execution times for different algorithms: Frechet-MLE")

##Plot of absolute times for talk----
###create plots----
####-factors----
benchTibMleN <- fullBenchTibN
benchTibMleN$expr <- benchTibMleN$expr %>% factor(
  levels = c("k2", "k3", "sliding", "k2No", "k3No","disjoint"), 
  labels = c("cb(2)", "cb(3)", "sb", "k2No", "k3No", "db")
)
benchTibMleN$B <- benchTibMleN$B %>% factor(
  levels = c(250, 500, 750, 1000), 
  labels = c("B = 250", "B = 500", "B = 750", "B = 1000")
)
####plotting----

absTimePlotN <-
  benchTibMleN  %>% 
  filter(!expr %in% c("k2No", "k3No")) %>% 
  ggplot(aes(x = m, y = timeMed/1000, col = expr))+ geom_line()+
  facet_wrap(~B, ncol = 4)+
  scale_colour_colorblind()+
  labs(x = "Effective sample size",
       title = "Absolute runtimes for Fréchet MLE", 
       y = "Time in seconds")+
  labsPlot+
  scale_y_continuous(limits = c(0,3.1), breaks = c(0,1,2,3))+
  themePlot+
  theme( 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x  = element_blank()
  )+ 
  labs(title = "Runtimes for Fréchet MLE")
absTimePlotN

#now see end of file for arranging!

##Plot of relative times----
###create plots----
####-factors----
tmpDbTibMleN <- fullBenchTibN %>% filter(expr == "disjoint") %>% rename(timeDB = "timeMed") %>% 
  select(-expr)
benchTibRelMleN <- left_join(fullBenchTibN, tmpDbTibMleN) %>% filter(expr !="disjoint") %>% 
  mutate(relTime = timeMed/timeDB)
benchTibRelMleN$expr <- benchTibRelMleN$expr %>% factor(
  levels = c("disjoint", "k2", "k3", "sliding", "k2No", "k3No"), 
  labels = c("db", "cb(2)", "cb(3)", "sb", "k2No", "k3No")
)
benchTibRelMleN$B <- benchTibRelMleN$B %>% factor(
  levels = c(250, 500, 750, 1000), 
  labels = c("B = 250", "B = 500", "B = 750", "B = 1000")
)
####plotting----
relTimePlotN <-
  benchTibRelMleN  %>% 
  filter(!expr %in% c("k2No", "k3No")) %>% 
  ggplot(aes(x = m, y = relTime, col = expr))+ geom_line()+
  scale_color_colorblind()+
  geom_hline(yintercept = 1, col = "black", linetype = "dashed")+
  facet_wrap(~B, ncol = 4)+
  labs(x = "Effective sample size",
       title = "Relative runtimes for Fréchet MLE", 
       y = "Relative time")+
  labsPlot+
  themePlot+
  theme( 
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
relTimePlotN
if(F){
ggsave(relTimePlotN, path = here("results/"), 
       filename = "plotrelTimePlotMleNFix.pdf", device = "pdf",
       width = 9, height = 3, 
)
}

###plots for talk----
relTimePlotTalkN <- relTimePlotN+
  theme(strip.text.x = element_blank(),
        plot.title = element_blank())
  
relTimePlotTalkN

talkPlot <- ggarrange( absTimePlotN, relTimePlotTalkN, align = "v", ncol = 1,
                       common.legend = T, legend = "bottom",
                       heights = c(1,1)
)
talkPlot
if(F){
ggsave(plot = talkPlot, path = here("results/"), 
       filename = "timesPlot.pdf", device = "pdf",
       width = 12, height = 6,
       )
}





