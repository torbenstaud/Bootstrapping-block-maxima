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
#Hardcoded variables
gammaVec <- c(-0.2, -0.1, 0, 0.1, 0.2)
gammaVecPub <- c(-0.2, 0, 0.2)
betaVec <- c(0, 0.5) #iid and ARMAX(0.5)
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







#bootstrap performance results----
#load array with CIs 
load(here(dataPath, "full_array_cis_ts.dat"))


#CI widths----
array_avg_ci_widths <- 
  apply(array_cis_full, c(1,2,3,4,5,6), function(xx) xx[3] - xx[2]) %>% 
  apply(c(2,3,4,5,6), mean)

#Transform to readable tibble
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
tibble_avg_ci_widths %>% 
  filter(parameter == plotParam) %>% 
  ggplot(
    aes(x = m, y = avgwidth, col = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "free", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  scale_color_manual(values = ownPalette)+
  themePlot+
  labs(title = paste0("Bootstrap confidence intervals for ", plotParam),
       y = "Average width", x = "Effective sample size m", col = "Method")

tibble_rel_avg_ci_widths <-
  tibble_avg_ci_widths %>% pivot_wider(names_from = method, values_from = avgwidth) %>% 
  mutate(sb = sb/db, `cb(2)` = `cb(2)`/db, db = 1) %>% 
  pivot_longer(cols = c(sb, `cb(2)`, db), names_to = "method") %>% 
  rename(rel_avg_width = value)

plotParam <- parameter_names[5]

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
ciRlWidthRelPlotJrssb

ciMeanWidthRelPlotJrssb <- 
  tibble_rel_avg_ci_widths %>% 
  filter(parameter == "mean", gamma %in% gammaVecPub, method !="sb") %>% 
  ggplot(aes(x = m, y = rel_avg_width, col = method, linetype = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "fixed", 
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
ciMeanWidthRelPlotJrssb
#now empirical coverage ----
##load truth
load(here(dataPath, "truth_array.dat"))


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

plotParam <- "Mean"

tibble_ci_covs$method <- factor(tibble_ci_covs$method,
                                levels = c("sb","db", "cb(2)"))
ciMeanCovPlotJrssb <-
  tibble_ci_covs%>% 
  filter(parameter == "mean", gamma %in% gammaVecPub, method !="sb") %>% 
  ggplot(
    aes(x = m, y = emp_cov, col = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "fixed", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  geom_abline(intercept = 0.95, slope = 0, col = "black", linetype = "longdash")+
  scale_color_manual(values = ownPalette, labels = c("db" ,"sb-cb"))+
  scale_y_continuous(limits = c(0.88, 1), breaks = c(0.9, 0.95, 1))+
  themePlot+
  labs(title = paste0("Bootstrap confidence intervals for ", plotParam),
       y = "Empirical coverage", x = "Effective sample size m", col = "Confidence Interval:")+
  theme(
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_blank()
  )
ciMeanCovPlotJrssb


combMeanCiPlotJrssb <- 
  ggarrange(ciMeanCovPlotJrssb, ciMeanWidthRelPlotJrssb, ncol = 2, align = "h", 
            common.legend = T, legend = "bottom",
            heights = c(1, 1.15))
combMeanCiPlotJrssb
if(T){
  ggsave(combMeanCiPlotJrssb, path = plotPath, 
         filename = "Figure F.4.pdf", device = "pdf",
         width = 13, height = 7
  )
}





#bootstrap variance estimator performance
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
 

meanBstVarPlot <- 
  tib_joined %>% 
  filter(
         parameter == "mean", gamma %in% c(-0.2, 0, 0.2)) %>% 
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
meanBstVarPlot  


if(T) {
  ggsave(meanBstVarPlot, path = plotPath, 
         filename = "Figure F.3.pdf", device = "pdf",
         width = 13, height = 7
  )
}


##size matched bootstrap confidence intervals
#load data: tib_full_simaCi_chars
###cutoffgam = 30 means no cutoff, 0.5 also available
cutoffGam <- 30
nameLoad <- sprintf("tib_full_cutoff%.1f_multSimaCiChars", cutoffGam)
load(file = file.path(dataPath, nameLoad))

tib_full_simaCi_chars$method <- factor(tib_full_simaCi_chars$method, 
                                       levels = c("sb", "cb(2)", "db"),
                                       labels = c("sb","sb-cb", "db"))




tib_full_simaCi_charsRel <- 
  tib_full_simaCi_chars %>% filter(char == "avg_width") %>% 
  pivot_wider(names_from = "method") %>% 
  mutate(sb = sb/db, `sb-cb`= `sb-cb`/db, db = 1) %>% 
  pivot_longer(cols = c(sb, db, `sb-cb`), names_to = "method") %>% 
  mutate(char = "rel_width") %>% 
  bind_rows(tib_full_simaCi_chars)



#factor correction included
load(file = file.path(dataPath, sprintf("coeff_rl%.1f", cutoffGam)))


##MEAN




##now do empirical size correction based on shape estimates
#RL check "factor_correcion.R"
#MEAN

get_corr_mean <- function(gamma, m, nometh = T){
  res_fac <- coeff_mean[1] + gamma*coeff_mean[2] + m*coeff_mean[3]
  if(nometh == F){
    res_fac <- res_fac + coeff_mean[4]
  }
  return(res_fac)
}
##gamma + m+ method+ method*m+method+gamma*m
get_corr_mean_alles <- function(gamma, m, method){
  #cb fall
  res_fac <- coeff_mean_alles[1] + gamma*coeff_mean_alles[2] + m*coeff_mean_alles[3] + gamma*m*coeff_mean_alles[6]
  #db fall
  res_fac <- ifelse(method == "db", res_fac + coeff_mean_alles[4] + m*coeff_mean_alles[5], res_fac)
  return(res_fac)
}

###load array with estimates inside
tib_shape_ests <- 
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
  filter(parameter %in% c("shape"), char == "est") %>% 
  select(-c(char, parameter)) %>% 
  rename(shape_est = value)

tib_cisMeanCorr_full <- 
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
  filter(parameter %in% c("mean")) %>% 
  pivot_wider(names_from = char) %>%
  full_join(tib_shape_ests) %>% 
  mutate(
    width = get_corr_mean(shape_est, m) * (upper - lower)) %>% 
  full_join(tibble_truth) %>%  
  mutate(
    lower = est + get_corr_mean(shape_est, m) *lower, 
    upper = est + get_corr_mean(shape_est, m) *upper, 
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


pltMeanNoMethFactorCorrectedEmpCov <- 
  tib_cisMeanCorr_full %>% 
  filter(parameter == parameter_names[5], gamma %in% gammaVecPub, method !="sb",
         char == "emp_cov") %>% 
  ggplot(
    aes(x = m, y = value, col = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "fixed", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  geom_abline(intercept = 0.95, slope = 0, col = "black", linetype = "longdash")+
  scale_color_manual(values = ownPalette, labels = c("db" ,"sb-cb"))+
  scale_y_continuous(limits = c(0.9, 1), breaks = c(0.9, 0.95, 1))+
  themePlot+
  labs(title = sprintf("Data adaptive size correction without method for mean"),
       y = "Empirical coverage", x = "Effective sample size m", col = "Confidence Interval:")+
  theme(
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_blank()
  )
pltMeanNoMethFactorCorrectedEmpCov 

pltMeanNoMethFactorCorrectedRelCi <-
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
  filter(parameter %in% c("mean")) %>% 
  pivot_wider(names_from = char) %>%
  full_join(tib_shape_ests) %>% 
  mutate(
    width = get_corr_mean(shape_est, m) * (upper - lower)) %>% 
  full_join(tibble_truth) %>%  
  mutate(
    lower = est + get_corr_mean(shape_est, m) *lower, 
    upper = est + get_corr_mean(shape_est, m) *upper, 
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


pltMeanNoMethFactorCorrectedRelCi <- 
  tib_cisMeanCorr_full %>% 
  filter(parameter == "mean", gamma %in% gammaVecPub, method !="sb",
         char == "relwidth") %>% 
  ggplot(
    aes(x = m, y = value, col = method))+
  facet_grid(rows = vars(gamma), cols = vars(beta), scales = "fixed", 
             labeller = labeller(gamma = facet_labels_gamma, 
                                 beta = facet_labels_beta))+
  geom_line(linewidth = 1)+
  scale_color_manual(values = ownPalette, labels = c("db" ,"sb-cb"))+
  themePlot+
  labs(title = sprintf("Data adaptive size correction without method for mean "),
       y = "Relative CI width", x = "Effective sample size m", col = "Confidence Interval:")+
  theme(
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_blank()
  )
pltMeanNoMethFactorCorrectedRelCi
combMeanCiPlotFaCoJrssb <- 
  ggarrange(pltMeanNoMethFactorCorrectedEmpCov, pltMeanNoMethFactorCorrectedRelCi, ncol = 2, align = "h",common.legend = T, legend = "bottom")
combMeanCiPlotFaCoJrssb

if(T){ 
  ggsave(combMeanCiPlotFaCoJrssb, path = plotPath, 
         filename = sprintf("Figure F.5.pdf"), device = "pdf",
         width = 13, height = 7
  )
}



