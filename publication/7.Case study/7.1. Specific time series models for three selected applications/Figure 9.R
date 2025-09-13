library(here)
library(tidyverse)
library(ggpubr)
textSize <- 20
source(here("backend/0plot_theme.R"))
section <- "7.Case study/7.1. Specific time series models for three selected applications"

dataPath <- here(file.path("backend", "data", "7", "7.1"))

parameter_names <- 
  c("location", "scale", "shape", "return_level(100,365)", "mean")



#Plotting of estimator performance----
facet_labels_model <- as_labeller(
  c(
    "log returns" = "Log Returns", 
    "precipitation" = "Precipitation",
    "temperature" = "Temperature"
  )
)

facet_labels_chars <- as_labeller(
  c(
    "emp_cov" = "Empirical coverage", 
    "rel_avg_width" = "Relative CI width",
    "rel_mse" = "Relative MSE",
    "bias_par" = "Bias^2/MSE",
    "mse" = "MSE"
  )
)
load(here(dataPath, "tibble_ests_prcp_365"))
load(here(dataPath, "tibble_ests_wxgen_92"))
load(here(dataPath, "tibble_ests_garch_250"))

tibble_ests_allmodels <-
  bind_rows(
    tibble_ests_prcp_365 %>% mutate(model = "precipitation", bias_sq = bias**2),
    tibble_ests_garch_250 %>% mutate(model = "log returns", bias_sq = bias**2),
    tibble_ests_wxgen_92 %>% mutate(model = "temperature", bias_sq = bias**2)
  ) %>% 
  pivot_longer(cols = c(bias, variance, mse, truth, emp_mean, bias_par, bias_sq),
               names_to = "char")
tibble_rel_mse_allmodels <- 
  tibble_ests_allmodels %>% filter(char == "mse") %>% 
  pivot_wider(names_from = method) %>% 
  mutate(`cb(2)` = `cb(2)`/db, sb = sb/db, db = 1, char = "rel_mse") %>% 
  pivot_longer(cols = c(`cb(2)`, sb, db), names_to = "method")
tibble_ests_allmodels <- full_join(tibble_ests_allmodels, tibble_rel_mse_allmodels)


tibble_ests_allmodels$char <- 
  factor(tibble_ests_allmodels$char, 
         levels = 
           c("mse","rel_mse", "bias_par", "bias_sq", "variance", "emp_mean", 
             "truth", "bias"))

tibble_ests_allmodels$method <- 
  factor(tibble_ests_allmodels$method, 
         levels = c( "cb(2)", "db", "sb"),
         labels = c("cb","db", "sb" ))



getPlotAppl <- function(model = "log returns", noyfac = F, nox = F){
  if(model == "log returns"){
    tib <- tibble_ests_allmodels %>% 
      filter(char %in% c("rel_mse", "mse"), 
             parameter == parameter_names[4], model == "log returns")
  } else if(model == "precipitation"){
    tib <- tibble_ests_allmodels %>% 
      filter(char %in% c("rel_mse", "mse"), 
             parameter == parameter_names[4], model == "precipitation")
  }else if(model == "temperature"){
    tib <- tibble_ests_allmodels %>% 
      filter(char %in% c("rel_mse", "mse"), 
             parameter == parameter_names[4], model == "temperature")
  }
  resPlt <- ggplot(tib, aes(x = m, y = value, col = method)) +
    facet_grid(rows = vars(char), cols = vars(model), scales = "free",
               labeller = labeller(model = facet_labels_model, 
                                   char = facet_labels_chars)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = ownPalette) +
    labs(x = "Effective sample size m",
         col = "Estimator:", linetype = "Confidence Interval")+
    themePlot +
    theme(axis.title.y = element_blank(),
          legend.position = "bottom")+
    theme(plot.margin = margin(t = 5, r = 0.5, b = 5, l = 0.5))
  if(noyfac == T){
    resPlt <- resPlt + 
      theme(strip.text.y = element_blank())
  }
  if(nox == T){
    resPlt <- resPlt + 
      xlab(" ")
  }
  return(resPlt)
}

pltApplEstPerf <- 
  ggarrange(
    getPlotAppl(model = "log returns", noyfac = T, nox = T),
    getPlotAppl(model = "precipitation", noyfac = T, nox = F),
    getPlotAppl(model = "temperature", noyfac =F, nox = T),
    nrow = 1, legend = "bottom", common.legend = T,
    align = "hv",
    widths = c(1, 1, 1)       # gleiche Breite
  )  
pltApplEstPerf
if(T){
  ggsave(pltApplEstPerf, 
         filename = "Figure 9.pdf",
         device = "pdf", path = file.path(here(section), "plots"),
         width = 13, height = 7)
}