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


textSize <- 20
source(here("backend/0plot_theme.R"))
plotPath <- file.path("Supplement", "plots")

dataPath <- here(file.path("backend", "data", "Supplement"))

# ==== SECTION 3: MAIN ANALYSIS ================================================
# Description: The core analytical processes, including data loading, cleaning,
#              processing, analysis, and visualization.
#______________________________________________________________________________


## ==== SECTION 3.1: fixed r (block size) ======================================
# Description: The core analytical processes, including data loading, cleaning,
#              processing, analysis, and visualization.
#_______________________________________________________________________________

#load the processed data file; name: fullBenchTibR
load(file.path(dataPath, "benchTibMleRFix"))


labsPlot <- labs(
  x = "Effective sample size m",
  col = "Bootstrap")
###-factors  ========================================================
tmpDbTibMle <- fullBenchTibR %>% filter(expr == "disjoint") %>% rename(timeDB = "timeMed") %>% 
  select(-expr)
benchTibRelMle <- left_join(fullBenchTibR, tmpDbTibMle) %>% filter(expr !="disjoint") %>% 
  mutate(relTime = timeMed/timeDB)
benchTibRelMle$expr <- benchTibRelMle$expr %>% factor(
  levels = c("k2", "sliding", "k3",  "k2No", "k3No","disjoint"), 
  labels = c("cb(2)", "sb", "cb(3)", "k2No", "k3No", "db")
)


benchTibRelMle$B <- benchTibRelMle$B %>% factor(
  levels = c(250, 500, 750, 1000), 
  labels = c("B = 250", "B = 500", "B = 750", "B = 1000")
)
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
###plotting----
benchTibRelMle$expr <- 
  factor(benchTibRelMle$expr, 
         levels = c("db", "sb", "cb(2)", "cb(3)", "k2No", "k3No"),
         labels = c("db", "sb", "cb", "cb(3)", "k2No", "k3No"))
relTimePlot <-
  benchTibRelMle  %>% mutate(n = m*r) %>% 
  filter(!expr %in% c("k2No", "k3No", "cb(3)"), m >= 40) %>% 
  ggplot(aes(x = n, y = relTime, col = expr))+ 
  geom_line(linewidth = 1.1)+
  geom_hline(yintercept = 1, col = "#00BFC4", linetype = "solid", linewidth = 1.1)+
  facet_wrap(~B, ncol = 4)+
  labsPlot+
  labs(x = "Total sample size n",
       title = "Relative runtimes for Fréchet MLE", 
       y = "Relative Runtime", 
       col = "Bootstrap:")+
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
if(TRUE){
  ggsave(relTimePlot, path = plotPath, 
         filename = "Figure F.11.pdf", device = "pdf",
         width = 13, height = 3.5, 
  )
}





