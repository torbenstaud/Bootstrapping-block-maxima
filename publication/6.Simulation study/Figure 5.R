library(tidyverse)
library(here)
library(Rcpp)
library(ggpubr)
library(tidyr)
library(dplyr)
section <- "6.Simulation study"

textSize <- 20
source(here("backend/0plot_theme.R"))
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




dataPath <- here(file.path("backend", "data", "6"))


#nun bootstrap Varianzschätzer Performance plotten
load(file.path(dataPath, "tib_full_bst_vars")) #tib_full_bst_vars
load(file.path(dataPath, "tib_full_estVars_ts")) #tib_full_true_estVar
tib_full_bst_vars <- tib_full_bst_vars %>% 
  mutate(target_var = ifelse(method == "cb(2)", "sb", method))
tib_full_true_estVar <- tib_full_true_estVar %>%
  rename(target_var = method, true_est_var = value)

tib_joined <- full_join(tib_full_bst_vars , tib_full_true_estVar)


facet_label_targetvar <- as_labeller(
  c(
    "sb" = "sigma[sb]^{2}", 
    "db" = "sigma[db]^{2}"
  ),
  label_parsed
)
labsPlot <- labs(
) 
tib_joined$target_var <- factor(tib_joined$target_var,
                                levels = c("sb", "db"))

rlBstVarPlot <- 
  tib_joined %>% 
  filter(method != "asb",
         parameter == "return_level(100,365)", gamma %in% c(-0.2, 0, 0.2)) %>% 
  ggplot(aes(x = m, y = value, col = method))+
  facet_grid(gamma ~ target_var + beta , scales = "free_y",
             labeller = labeller( target_var = facet_label_targetvar, 
                                  gamma = facet_labels_gamma, 
                                  beta = facet_labels_beta))+
  geom_line(linewidth = 1.1)+
  geom_line(aes(x = m, y = true_est_var), 
            col = "black", linetype = "longdash",
            linewidth = 1.1)+
  scale_color_manual(values = ownPalette, labels = c("cb", "db", "naive sb"))+
  labs( y = paste0("Variance"),
        x = "Effective sample size m",
        col = "Bootstrap:")+
  themePlot+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title = element_blank())
rlBstVarPlot  


if(T) {
  ggsave(rlBstVarPlot, path = file.path(here(section), "plots"), 
         filename = "Figure 5.pdf", device = "pdf",
         width = 13, height = 7 #urspr 10, 12
  )
}
