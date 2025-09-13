library(here)
library(tidyverse)
textSize <- 20
source(here("backend/0plot_theme.R"))
section <- "7.Case study/7.1. Specific time series models for three selected applications"

dataPath <- here(file.path("backend", "data", "7", "7.1"))

parameter_names <- 
  c("location", "scale", "shape", "return_level(100,365)", "mean")


#plotting for bootstrap performance----


load(here(dataPath, "tibble_corrected_rel_cis_Garch"))
load(here(dataPath, "tibble_corrected_rel_cis_viet"))
load(here(dataPath, "tibble_corrected_rel_cis_wgen"))

tibble_full <- 
  bind_rows(tibble_corrected_rel_cis_Garch, 
            tibble_corrected_rel_cis_viet, 
            tibble_corrected_rel_cis_wgen)

facet_labels_chars <- as_labeller(
  c(
    "emp_cov" = "Empirical coverage", 
    "rel_avg_width" = "Relative CI width",
    "rel_mse" = "Relative MSE",
    "bias_par" = "Bias^2/MSE",
    "mse" = "MSE"
  )
)
facet_labels_model <- as_labeller(
  c(
    "log returns" = "Log Returns", 
    "precipitation" = "Precipitation",
    "temperature" = "Temperature"
  )
)
tibble_full$method <- factor(tibble_full$method, levels = c("sb", "db", "cb(2)"),
                             labels = c("sb", "db", "sb-cb"))
plotCommon <- 
  tibble_full %>% ggplot(aes(x = m, y = value, col = method)) +
  facet_grid(rows = vars(char), cols = vars(model), scales = "free",
             labeller = labeller(model = facet_labels_model, char = facet_labels_chars)) +
  geom_line(linewidth = 1.1) +
  
  # Fixiere y-Achse für "emp_cov" auf [0.9, 1]
  geom_blank(data = tibble_full %>% filter(char == "emp_cov") %>% mutate(value = 0.9)) + 
  geom_blank(data = tibble_full %>% filter(char == "emp_cov") %>% mutate(value = 1.0)) +
  
  # geom_abline für "emp_cov" nur in der ersten Zeile
  geom_abline(data = tibble_full %>% 
                filter(char == "emp_cov") %>% mutate(dummy_x = NA, dummy_y = NA), 
              aes(intercept = 0.95, slope = 0), color = "black", linetype = "longdash") + 
  
  scale_color_manual(values = ownPalette) +
  labs(x = "Effective sample size m",
       col = "Confidence Interval", linetype = "Confidence Interval") +
  themePlot +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")
plotCommon

if(T){
  ggsave(plotCommon, 
         filename = "Figure 10.pdf",
         device = "pdf", path = file.path(here(section), "plots"),
         width = 13, height = 7)
}


