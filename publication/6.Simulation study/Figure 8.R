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


parameter_names <- c("location", "scale", "shape", "return_level(100,365)", "mean", 
                     "return_level(75,365)", "return_level(50,365)")
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

#factor correction script
##load size matched data: tib_full_simaCi_chars
cutoffGam <- 30  #30(so no cutoff) #0.5 also available
nameLoad <- sprintf("tib_full_cutoff%.1f_multSimaCiChars", cutoffGam)
load(file = file.path(dataPath, nameLoad))

#derive coefficients for the factors for RL parameter
tib_cleaned_factors <- 
  tib_full_simaCi_chars %>% 
  filter(char == "factor", method != "sb", parameter==parameter_names[4]) %>% 
  select(-c(parameter, char))

mod_rl <- lm(value ~ gamma + beta +m + method, data = tib_cleaned_factors) 

coeff_rl <- (mod_rl$coeff)[c(1,2,4)]
tib_c_plot <- 
  bind_rows(
    tibble(m = seq(40, 100, length.out = 10**3), gamma = -0.2) %>% 
      mutate(c = m*coeff_rl[3] + coeff_rl[2]*gamma + coeff_rl[1]),
    tibble(m = seq(40, 100, length.out = 10**3), gamma = 0) %>% 
      mutate(c = m*coeff_rl[3] + coeff_rl[2]*gamma + coeff_rl[1]),
    tibble(m = seq(40, 100, length.out = 10**3), gamma = 0.2) %>% 
      mutate(c = m*coeff_rl[3] + coeff_rl[2]*gamma + coeff_rl[1])
  )

tib_c_plot$gamma <- factor(tib_c_plot$gamma)
plot_rl_factors <- 
  tib_c_plot %>% 
  ggplot(aes(x = m, y = c, linetype = gamma))+
  #facet_grid(cols = vars(gamma), 
  #           labeller = labeller(gamma = facet_labels_gamma ))+
  geom_line(linewidth = 1.1)+
  labs( y = "Correction Factor c",
        x = "Effective sample size m",
        linetype = expression(gamma))+
  themePlot+
  theme(plot.title = element_blank(),
        strip.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#derive coefficients for the factors for mean parameter
tib_cleaned_factors <- 
  tib_full_simaCi_chars %>% 
  filter(char == "factor", method != "sb", parameter==parameter_names[5]) %>% 
  select(-c(parameter, char))

mod_rl <- lm(value ~ gamma + m+ method, data = tib_cleaned_factors) 
#mod_rl %>% summary
coeff_mean <- (mod_rl$coeff)[c(1,2,3)]
coeff_mean_alles <- mod_rl$coeff

save(coeff_rl, coeff_mean, coeff_mean_alles, file = file.path(dataPath, sprintf("coeff_rl%.1f", cutoffGam)))
#not differentiating between (block maximum) method
mod_rl <- lm(value ~ gamma + beta +m , data = tib_cleaned_factors) 
mod_rl %>% summary
coeff_mean <- (mod_rl$coeff)[c(1,2,4)]
save(coeff_rl, coeff_mean, file = file.path(dataPath, sprintf("coeff_rl_nometh%.1f", cutoffGam)))


# factor correction for time-series model

get_corr <- function(gamma, m){
  coeff_rl[1] + gamma*coeff_rl[2] + m*coeff_rl[3]
}
#loading files
load(file.path(dataPath, "truth_array.dat"))
load(file.path(dataPath, "full_array_cis_ts.dat"))

#Transform to humanly interpretable tibble
parameter_names <- c("location", "scale", "shape", "return_level(100,365)", "mean", 
                     "return_level(75,365)", "return_level(50,365)")
method_names <- c("sb", "db", "cb(2)")

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
  select(gamma, beta, parameter, value) %>% rename(truth = value) %>% 
  filter( parameter == parameter_names[4])


tib_cisCorr_full <- 
  as_tibble(as.data.frame(as.table(array_cis_full))) %>%
  rename(N_index = Var1, gamma_index = Var2, beta_index = Var3, m_index = Var4, 
         parameter_index = Var5, method_index = Var6, char_index = Var7,
         value = Freq) %>%
  mutate(
    N_index = seq(1, N)[N_index],
    gamma = gammaVec[gamma_index],
    beta = betaVec[beta_index],
    m = mVec[m_index],
    parameter = parameter_names[parameter_index], 
    char = (c("est", "lower", "upper"))[char_index],
    method = method_names[method_index]
  ) %>%
  select(N_index, gamma, beta, m, parameter, char, method, value) %>% 
  filter(parameter == parameter_names[4]) %>% 
  pivot_wider(names_from = char) %>% 
  mutate(
    width = get_corr(gamma, m) * (upper - lower), 
    parameter = "rl"
  ) %>% 
  full_join(tibble_truth) %>%  
  mutate(
    lower = est + get_corr(gamma, m) *lower, 
    upper = est + get_corr(gamma, m) *upper, 
    cov = ifelse(lower <= truth & truth <= upper, 1, 0)
  ) %>% 
  as_tibble() %>%  
  group_by(gamma, beta, m, parameter, method) %>% 
  summarise(avgwidth = mean(width), emp_cov = mean(cov)) %>% 
  select(gamma, beta, m, parameter, method, avgwidth, emp_cov) %>% 
  filter(gamma %in% c(-0.2, 0, 0.2), method != "sb") %>% 
  pivot_longer(cols = c("avgwidth", "emp_cov"), names_to = "char") %>% 
  pivot_wider(names_from = method) %>% 
  mutate(`cb(2)` = ifelse(char == "avgwidth", `cb(2)`/db, `cb(2)`),
         db = ifelse(char == "avgwidth", 1, db)) %>% 
  pivot_longer(cols = c(`cb(2)`, db), names_to = "method") %>% 
  mutate(char = ifelse(char == "avgwidth", "relwidth", char))




plotParam <- parameter_names[4]

tib_cisCorr_full <- tib_cisCorr_full %>% 
  mutate(method = ifelse(method == "cb(2)", "sb-cb", method))
pltFaCoRelWidthTs <- tib_cisCorr_full %>% 
  filter(parameter == "rl", gamma %in% gammaVecPub, method !="sb", char == "emp_cov") %>% 
  ggplot(
    aes(x = m, y = value, col = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "free", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  geom_abline(intercept = 0.95, slope = 0, col = "black", linetype = "longdash")+
  scale_color_manual(values = ownPalette)+
  scale_y_continuous(limits = c(0.9, 1), breaks = c(0.9, 0.95, 1))+
  themePlot+
  labs(title = paste0("Bootstrap confidence intervals for ", plotParam),
       y = "Empirical coverage", x = "Effective sample size m", col = "Confidence Interval:")+
  theme(
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_blank()
  )

pltFaCoEmpCovTs <- tib_cisCorr_full %>% 
  filter(parameter == "rl", gamma %in% gammaVecPub, 
         method !="sb", char == "relwidth") %>% 
  ggplot(
    aes(x = m, y = value, col = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "free", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  scale_color_manual(values = ownPalette)+
  scale_linetype_manual(values = c("solid", "solid", "solid"))+
  themePlot+
  labs(title = paste0("Bootstrap confidence intervals for ", plotParam),
       y = "Relative average width", x = "Effective sample size m", col = "Confidence Interval:",
       linetype = "Confidence Interval:")+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



#combine plots
combRlFaCoCiPlotTsJrssb <- 
  ggarrange(pltFaCoRelWidthTs, pltFaCoEmpCovTs, ncol = 2, align = "h", 
            common.legend = T, legend = "bottom",
            heights = c(1, 1.15))
combRlFaCoCiPlotTsJrssb

if(T){
  ggsave(combRlFaCoCiPlotTsJrssb, path = file.path(here(section), "plots"), 
         filename = "Figure 8.pdf", device = "pdf",
         width = 13, height = 7
  )
}
