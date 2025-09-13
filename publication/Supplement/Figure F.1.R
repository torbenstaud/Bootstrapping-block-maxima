library(tidyverse)
library(here)
library(Rcpp)
library(ggpubr)
library(tidyr)
library(dplyr)

textSize <- 20
source(here("backend/0plot_theme.R"))
plotPath <- file.path("Supplement", "plots")

dataPath <- here(file.path("backend", "data", "Supplement"))
#Alle harcoded variablen definieren
gammaVec <- c(-0.2, -0.1, 0, 0.1, 0.2)
gammaVecPub <- c(-0.2, 0, 0.2)
betaVec <- c(0, 0.5) #iid und einmal ARMAX(0.5)
mVec <- seq(40,100, by = 10)
r <- 365
N <- 10**3
B <- 10**3
NSeedChunkSize <- 150
NSeedChunkN <- N/NSeedChunkSize
NSeedChunkVec <- seq(1, NSeedChunkN)
Time <- 100
level <- 0.05


#estimator performance results----
#load data
##load estimates
load(file.path(dataPath,"full_tib_ests"))
tib_ests <- tib_ests %>% 
  pivot_longer(cols = c(loc, scale, shape, rl, mean), names_to = "parameter")
## load true values
load(file.path(dataPath, "truth_array.dat"))

parameter_names <- c("loc", "scale", "shape", "rl", "mean")
tibble_truth <- 
  as_tibble(as.data.frame(as.table(array_truth))) %>% 
  rename(gamma_index = Var1, beta_index = Var2, parameter_index = Var3, 
         value = Freq)%>%
  mutate(
    gamma = gammaVec[gamma_index],
    beta = betaVec[beta_index],
    parameter = parameter_names[parameter_index],
  ) %>%
  select(gamma, beta, parameter, value) %>% rename(truth = value)

##join both estimates and true values
tib_full <- full_join(tib_ests, tibble_truth)

tib_chars <- 
  tib_full %>% mutate(diff = value - truth) %>% 
  group_by(gamma, beta, m, method, parameter) %>% 
  summarise(bias_sq = (mean(diff))**2, variance = var(value)) %>% 
  mutate(mse = bias_sq + variance) %>% 
  pivot_longer(cols = c(bias_sq, variance, bias_sq, mse), names_to = "char") 

facet_labels_gamma <- as_labeller(
  c(
    "-0.2" = "gamma == -0.2", 
    "-0.1" = "gamma == -0.1", 
    "0"    = "gamma == 0", 
    "0.1"  = "gamma == 0.1", 
    "0.2"  = "gamma == 0.2"
  ),
  label_parsed
)
facet_labels_beta <- 
  as_labeller(
    c("0" = "i.i.d.",
      "0.5" = "ARMAX 0.5"))



get_mse_relmse_arranged_plot <- function(plotParam, scale = 1){
  if(scale != 1){
    y_text <- paste0("MSE *", scale)
  }else{
    y_text <- "MSE"
  }
  pltRlEsts <- 
    tib_chars %>% 
    filter( char == "mse", parameter == plotParam, gamma %in% gammaVecPub) %>% 
    ggplot(
      aes(x = m, y = value*scale, col = method))+
    facet_grid(rows = vars(gamma), cols = vars(beta), scales = "free", 
               labeller = labeller(gamma = facet_labels_gamma, 
                                   beta = facet_labels_beta))+
    geom_line(linewidth = 1)+
    scale_color_manual(values = ownPalette)+
    labs( y = y_text,
          x = "Effective sample size m",
          col = "Estimator:")+
    themePlot+
    theme(plot.title = element_blank(),
          strip.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  pltRlEsts
  
  
  tib_chars_rel <- 
    tib_chars %>% pivot_wider(names_from = method) %>% 
    mutate(sb = sb/db, cb = cb/db, db = 1) %>% 
    pivot_longer(cols = c(sb, db, cb), names_to = "method")
  
  pltRlRelEsts <-
    tib_chars_rel %>% 
    filter( char == "mse", parameter == plotParam, gamma %in% gammaVecPub) %>% 
    ggplot(
      aes(x = m, y = value, col = method, linetype = method))+
    facet_grid(rows = vars(gamma), cols = vars(beta),  
               labeller = labeller(gamma = facet_labels_gamma, 
                                   beta = facet_labels_beta))+
    geom_line(linewidth = 1.1)+
    scale_color_manual(values = ownPalette)+
    scale_linetype_manual(values = c("solid", "solid", "solid"))+
    labs( y = "Relative MSE",
          x = "Effective sample size m",
          #title = paste0(plotParam," für wxgenr estimators; r = 365"), 
          col = "Estimator:", linetype = "Estimator:")+
    themePlot+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_blank())
  pltRlRelEsts
  
  plotAbsRelMseJrssb <- 
    ggarrange(pltRlEsts, pltRlRelEsts, align = "h", ncol = 2,
              common.legend = T, legend = "bottom")
  return(plotAbsRelMseJrssb)
}
plotRl <-
  get_mse_relmse_arranged_plot("rl", 1)
plotRl
if(T){
  ggsave(plotRl, filename = "Figure F.1.pdf",
         device = "pdf", path = here(plotPath),
         width = 13, height = 7)
}
plotMean <- 
  get_mse_relmse_arranged_plot(parameter_names[5], 100)
plotMean
if(FALSE){
  ggsave(plotMean, filename = "plotAbsRelMseMeanJrssb.pdf",
         device = "pdf", path = here("results/"),
         width = 13, height = 7)
}