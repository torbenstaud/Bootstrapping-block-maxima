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

#bootstrap performance results----
#array mit CIs laden
load(file.path(dataPath,"full_array_cis_ts.dat"))


#CI widths----
array_avg_ci_widths <- 
  apply(array_cis_full, c(1,2,3,4,5,6), function(xx) xx[3] - xx[2]) %>% 
  apply(c(2,3,4,5,6), mean)

#Transformiere in lesbares tibble
parameter_names <- c("location", "scale", "shape", "return_level(100,365)", "mean", 
                     "return_level(75,365)", "return_level(50,365)")
method_names <- c("sb", "db", "cb(2)")

tibble_avg_ci_widths <- as_tibble(as.data.frame(as.table(array_avg_ci_widths))) %>%
  rename(gamma_index = Var1, beta_index = Var2, m_index = Var3, 
         parameter_index = Var4, method_index = Var5, 
         avgwidth = Freq) %>%
  mutate(
    gamma = gammaVec[gamma_index],
    beta = betaVec[beta_index],
    m = mVec[m_index],
    parameter = parameter_names[parameter_index], 
    method = method_names[method_index]
  ) %>%
  select(gamma, beta, m, parameter, method, avgwidth) %>% 
  filter(gamma %in% c(-0.2, 0, 0.2))

plotParam <- parameter_names[5]

tibble_rel_avg_ci_widths <-
  tibble_avg_ci_widths %>% pivot_wider(names_from = method, values_from = avgwidth) %>% 
  mutate(sb = sb/db, `cb(2)` = `cb(2)`/db, db = 1) %>% 
  pivot_longer(cols = c(sb, `cb(2)`, db), names_to = "method") %>% 
  rename(rel_avg_width = value)


tibble_rel_avg_ci_widths$method <- factor(tibble_rel_avg_ci_widths$method,
                                          levels = c("sb","db", "cb(2)"))
ciRlWidthRelPlotJrssb <- 
  tibble_rel_avg_ci_widths %>% 
  filter(parameter == parameter_names[4], gamma %in% gammaVecPub, method !="sb") %>% 
  ggplot(aes(x = m, y = rel_avg_width, col = method, linetype = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "fixed", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  scale_color_manual(values = ownPalette)+
  scale_linetype_manual(values = c("solid", "solid", "solid"))+
  #scale_y_continuous(limits = c(0.75,1))+
  themePlot+
  labs(title = paste0("Bootstrap confidence intervals for ", plotParam),
       y = "Relative average width", x = "Effective sample size m", col = "Confidence Interval:",
       linetype = "Confidence Interval:")+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#nun empirical coverage ----
##zunächst wahrheit laden
load(file.path(dataPath,"truth_array.dat"))

#name: array_cis_full; dim [N, 5, 2, 7, 5, 3, 3]
array_cis <- array(dim = c(N, 5, 2, 7, 5, 3,2))
array_cis_full[,,,,,3,1] <- array_cis_full[,,,,,1,1]#set cb senter to sb point estimate
array_cis[,,,,,,1] <-   array_cis_full[,,,,, ,1] + array_cis_full[,,,,, , 2]
array_cis[,,,,,,2] <-   array_cis_full[,,,,, ,1] + array_cis_full[,,,,, , 3]

expand_truth <- array(dim = c(N, 5, 2, 7, 5, 3))
for(indN in seq(1,N)){
  for(indGamma in seq_along(gammaVec)){
    for(indBeta in seq_along(betaVec)){
      for(indM in seq_along(mVec)){
        for(indParam in seq(1,5)){
          for(indMeth in seq(1,3)){
            expand_truth[indN, indGamma, indBeta, indM, indParam, indMeth] <-
              array_truth[indGamma, indBeta, indParam]
          }
        }
      }
    }
  }
}
calc_empcov <- function(array_ci_truth){
  #first two: lower and upper CI interval end points
  #third: truth
  #array_ci dim(N, 3), truth is a number
  #array_ci <- 
  #  array_ciNums[,1] + array_ciNums[,c(2,3)]
  empCov <- 
    apply(array_ci_truth, c(1), function(xx) c(xx[1] <= xx[3] & xx[3] <= xx[2], 
                                               xx[1] > xx[3], xx[2] < xx[3])) %>% 
    aperm(c(2,1))
  return(apply(empCov, c(2), mean))
  
}
array_cis_truth <- array(dim = c(N, 5, 2, 7, 5, 3,3))
array_cis_truth[, , , , , ,c(1,2)] <- array_cis
array_cis_truth[, , , , , ,3] <- expand_truth


array_ci_covs <- array_cis_truth %>% apply(
  c(2,3,4,5,6), calc_empcov
) %>% aperm(c(2,3,4,5,6,1))
array_ci_covs[5,2,7,5,,1]


tibble_ci_covs <- 
  as_tibble(as.data.frame(as.table(array_ci_covs[,,,,,1]))) %>%
  rename(gamma_index = Var1, beta_index = Var2, m_index = Var3, 
         parameter_index = Var4, method_index = Var5, 
         emp_cov = Freq) %>%
  mutate(
    gamma = gammaVec[gamma_index],
    beta = betaVec[beta_index],
    m = mVec[m_index],
    parameter = parameter_names[parameter_index], 
    method = method_names[method_index]
  ) %>%
  select(gamma, beta, m, parameter, method, emp_cov)

plotParam <- "return_level(100,365)"

tibble_ci_covs$method <- factor(tibble_ci_covs$method,
                                levels = c("sb","db", "cb(2)"))
ciRlCovPlotJrssb <-
  tibble_ci_covs%>% 
  filter(parameter == parameter_names[4], gamma %in% gammaVecPub, method !="sb") %>% 
  ggplot(
    aes(x = m, y = emp_cov, col = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "fixed", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  geom_abline(intercept = 0.95, slope = 0, col = "black", linetype = "longdash")+
  scale_color_manual(values = ownPalette, labels = c("db" ,"sb-cb"))+
  scale_y_continuous(limits = c(0.7, 1))+
  themePlot+
  labs(title = paste0("Bootstrap confidence intervals for ", plotParam),
       y = "Empirical coverage", x = "Effective sample size m", col = "Confidence Interval:")+
  theme(
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_blank()
  )

combRlCiPlotJrssb <- 
  ggarrange(ciRlCovPlotJrssb, ciRlWidthRelPlotJrssb, ncol = 2, align = "h", 
            common.legend = T, legend = "bottom",
            heights = c(1, 1.15))

if(T){
  ggsave(combRlCiPlotJrssb, path = file.path(here(section), "plots"), 
         filename = "Figure 6.pdf", device = "pdf", 
         width = 13, height = 7
  )
}
