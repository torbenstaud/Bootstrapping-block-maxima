library(tidyverse)
library(here)
library(Rcpp)
library(ggpubr)
library(tidyr)
library(dplyr)
section <- "6.Simulation study"

textSize <- 20
source(here("backend/0plot_theme.R"))
source(here("backend", "data", "6", "labels.R"))
dataPath <- here(file.path("backend", "data", "6"))


gammaVec <- c(-0.2, -0.1, 0, 0.1, 0.2)
gammaVecPub <- c(-0.2, 0, 0.2)
betaVec <- c(0, 0.5) 
mVec <- seq(40,100, by = 10)
r <- 365
N <- 10**3
B <- 10**3
NSeedChunkSize <- 150
NSeedChunkN <- N/NSeedChunkSize
NSeedChunkVec <- seq(1, NSeedChunkN)
Time <- 100
level <- 0.05


cutoffGam <- 30
nameLoad <- sprintf("tib_full_cutoff%.1f_multSimaCiChars", cutoffGam)
load(file = file.path(dataPath, nameLoad))

tib_full_simaCi_chars$method <- factor(tib_full_simaCi_chars$method, 
                                       levels = c("sb", "cb(2)", "db"),
                                       labels = c("sb","sb-cb", "db"))






#plot correction factors for the RL(100) bootstrapping
plotParam <- "return_level(100,365)"
tib_full_simaCi_charsRel <- 
  tib_full_simaCi_chars %>% filter(char == "avg_width") %>% 
  pivot_wider(names_from = "method") %>% 
  mutate(sb = sb/db, `sb-cb`= `sb-cb`/db, db = 1) %>% 
  pivot_longer(cols = c(sb, db, `sb-cb`), names_to = "method") %>% 
  mutate(char = "rel_width") %>% 
  bind_rows(tib_full_simaCi_chars)


plotSimaQuantFac <- 
  tib_full_simaCi_charsRel  %>% 
  filter(char == "factor", parameter == plotParam, method != "sb",
         gamma %in%gammaVecPub) %>% 
  ggplot(aes(x = m, y = value, col = method))+
  facet_grid(gamma~beta, scales = "fixed",
             labeller = labeller( target_var = facet_label_targetvar, 
                                  gamma = facet_labels_gamma, 
                                  beta = facet_labels_beta))+
  geom_line(linewidth = 1.1)+
  scale_color_manual(values = ownPalette)+
  labs( y = "Correction Factor c for",
        x = "Effective sample size m",
        col = "CI Anchor",
        title = sprintf("Correction factor c for %s with cutoffgam = %.1f", plotParam, cutoffGam))+
  themePlot+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right"
  )


#plot relative average widths of RL cis
parameter_names <- c("location", "scale", "shape", "return_level(100,365)", "mean", 
                     "return_level(75,365)", "return_level(50,365)")
method_names <- c("sb", "db", "cb(2)")
plotChar <- "rel_width"
cap <- paste0(plotChar)
tib_full_simaCi_charsRel <- tib_full_simaCi_charsRel %>%  
  filter(method != "sb",
         gamma %in%gammaVecPub)
plotParam <- parameter_names[4]
plotSimaRlCov <- 
  tib_full_simaCi_charsRel  %>% 
  filter(char == "rel_width", parameter == plotParam, method != "sb",
         gamma %in% gammaVecPub) %>% 
  ggplot(aes(x = m, y = value, col = method))+
  facet_grid(gamma~beta, scales = "fixed",
             labeller = labeller( target_var = facet_label_targetvar, 
                                  gamma = facet_labels_gamma, 
                                  beta = facet_labels_beta))+
  geom_line(linewidth = 1.1)+
  scale_color_manual(values = ownPalette)+
  labs( y = "Relative CI width",
        x = "Effective sample size m",
        col = "Confidence Interval:",
        #title = sprintf("Relative CI width for %s with cutoffgam = %.1f", plotParam, cutoffGam)
        )+
  scale_y_continuous(limits = c(0.55, 1.05))+
  themePlot+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right"
  )




#factor correction included
load(file = file.path(dataPath, sprintf("coeff_rl%.1f", cutoffGam)))

facet_lines <- tib_full_simaCi_charsRel %>% filter(gamma %in% gammaVecPub) %>% 
  distinct(gamma, beta) %>%  
  mutate(
    intercept = coeff_rl[1] + coeff_rl[2] * gamma,  
    slope = coeff_rl[3]  
  )



plotParam <- parameter_names[4]
plotSimaQuantCovFct <- 
  tib_full_simaCi_charsRel %>%
  filter(char == "factor", parameter == "return_level(100,365)", method != "sb",
         gamma %in% gammaVecPub) %>%
  ggplot(aes(x = m, y = value, col = method)) +
  
  
  facet_grid(gamma ~ beta, scales = "fixed",
             labeller = labeller(target_var = facet_label_targetvar, 
                                 gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta)) +
  
  
  geom_line(linewidth = 1.1) +
  
  
  geom_abline(data = facet_lines, aes(intercept = intercept, slope = slope), 
              color = "black", linetype = "longdash") +
  
  
  scale_color_manual(values = ownPalette) +
  
  
  labs(y = "Correction Factor c",
       x = "Effective sample size m",
       col = "Confidence Interval",
       #title  = sprintf("Correction factor c for %s with cutoffgam = %.1f", plotParam, cutoffGam)
       ) +
  

  themePlot +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right")
plotSimaQuantCovFct



plotSimaRelCiFac <- 
  ggarrange(plotSimaRlCov, plotSimaQuantCovFct, legend = "bottom", 
            common.legend = T)

if(T){
  ggsave(plotSimaRelCiFac, path = file.path(here(section), "plots"), 
         filename = sprintf("Figure 7.pdf"), device = "pdf",
         width = 13, height = 7
  )
}
