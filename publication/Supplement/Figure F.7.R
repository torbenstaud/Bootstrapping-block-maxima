library(here)
library(tidyverse)
textSize <- 20
source(here("backend/0plot_theme.R"))
plotPath <- file.path("Supplement", "plots")

dataPath <- here(file.path("backend", "data", "Supplement", "F.2."))



#Hardcoded variables
copVec <- c("indep", "gauss", "gumbel", "tcop")
betaVec <- c(0, 0.5)
rhoVec <- c(0.2, 0.4) 
mVec <- seq(40,100, by = 10)
r <- 365
N <- 10**3
B <- 10**3
NSeedChunkSize <- 500
NSeedChunkN <- N/NSeedChunkSize
NSeedChunkVec <- seq(1, NSeedChunkN)
Time <- 100
level <- 0.05

#load array with CIs
load(here(dataPath, "full_array_cis_spearman.dat"))

#dim [N, 4, 2, 2, 7, 3, 3]; 1:Reps, 2:cop, 3:rho, 4:beta, 5:m, 6:method, 7:est value; estvalue - 0.975q; estvalue - 0.025q

#CI widths----
array_avg_ci_widths <- 
  apply(array_cis_full, c(1,2,3,4,5,6), function(xx) xx[3] - xx[2]) %>% 
  apply(c(2,3,4,5,6), mean)

#Transform to readable tibble
method_names <- c("sb", "db", "sb-cb")

tibble_avg_ci_widths_spearman <- as_tibble(as.data.frame(as.table(array_avg_ci_widths))) %>%
  rename(cop_index = Var1, rho_index = Var2, beta_index = Var3, m_index = Var4, 
         method_index = Var5, 
         avgwidth = Freq) %>%
  mutate(
    cop = copVec[cop_index],
    rho = rhoVec[rho_index],
    beta = betaVec[beta_index], 
    m = mVec[m_index],
    method = method_names[method_index]
  ) %>%
  select(cop, rho, beta, m, method, avgwidth)

facet_labels_cop <- as_labeller(
  c(
    "indep" = "Independence", 
    "gauss" = "Gauss", 
    "gumbel"    = "Gumbel-Hougaard", 
    "tcop"  = "t[4]"
  ),
  label_parsed
)

facet_labels_rho <- 
  as_labeller(
    c("0.2" = "rho == 0.2",
      "0.4" = "rho == 0.4"), 
    label_parsed)

tibble_avg_ci_widths_spearman %>% 
  filter(beta == 0.5, method != "sb") %>% 
  ggplot(
    aes(x = m, y = avgwidth, col = method))+
  facet_grid(rows = vars(rho), cols = vars(cop), 
             labeller = labeller(cop = facet_labels_cop, 
                                 rho = facet_labels_rho))+
  geom_line(linewidth=1.1)+
  scale_color_manual(values = ownPalette)+
  themePlot+
  labs(title = paste0("Bootstrap confidence intervals for Spearman's rho"),
       y = "Average width", x = "Effective sample size m", col = "Method")

tibble_rel_avg_ci_width_spearman <-
  tibble_avg_ci_widths_spearman %>% 
  pivot_wider(names_from = method, values_from = avgwidth) %>% 
  mutate(sb = sb/db, `sb-cb` = `sb-cb`/db, db = 1) %>% 
  pivot_longer(cols = c(sb, db, `sb-cb`), names_to = "method") %>% 
  rename(rel_avg_width = value) %>% 
  mutate(rho = ifelse(cop == "indep", 0, rho)) %>% 
  filter(method != "sb")


tibble_rel_avg_ci_width_spearman$rho <- 
  factor(tibble_rel_avg_ci_width_spearman$rho)

tibble_rel_avg_ci_width_spearman$cop <- 
  factor(tibble_rel_avg_ci_width_spearman$cop, 
         levels = 
           c("indep", "gauss", "tcop", "gumbel"))
tibble_rel_avg_ci_width_spearman <- 
  bind_cols(tibble_rel_avg_ci_width_spearman, tibble( char ="rel_avg_width")) %>% 
  rename(value = rel_avg_width)

tibble_rel_avg_ci_width_spearman %>% 
  filter(method != "sb", beta == 0) %>% 
  ggplot(aes(x = m, y = value, col = method, linetype = rho))+
  #scale_y_continuous(limits = c(0.9, 1.8))+
  facet_grid(cols = vars(cop), 
             labeller = labeller(cop = facet_labels_cop), scales = "free_y")+ 
  geom_line(linewidth=1.1)+
  geom_abline(slope = 0, intercept = 1, col = color_db, linewidth = 1.1)+
  scale_color_manual(values = ownPalette)+
  scale_linetype_manual(values = c("dotted", "solid", "longdash"))+
  themePlot+
  labs(title = "Bootstrap confidence intervals",
       y = "Relative average width", x = "Effective sample size m",
       col = "Confidence Intervals", linetype = expression(rho))+
  theme(legend.position = "bottom")

#Truths
tib_truthTemp <- tibble()
for(indCop in seq_along(copVec)){
  cop <- copVec[indCop]
  for(indRho in seq_along(rhoVec)){
    rho <- rhoVec[indRho]
    for(indBeta in seq_along(betaVec)){
      beta <- betaVec[indBeta]
      nameLoad <-
        paste0(
          "true_gen_cop", cop, "_rho", rho, "_beta", sprintf("%.1f",beta),".rds"
        )
      load(here(dataPath, nameLoad))
      for(indM in seq_along(mVec)){
        m <- mVec[indM]
        tib_truthTemp <- 
          bind_rows(tib_truthTemp,
                    bind_cols(tibble_res, tibble(m = m))
          )
        
      }
    }
  }
}

#Empirical Coverage
tib_truth <- 
  tib_truthTemp %>% 
  rename(rho = rho_innov, truth = rho_bm, cop = coptype) %>% 
  select(cop, beta, rho, m, truth) %>% 
  filter(!(cop == "indep" & rho == 0.4)) %>% 
  mutate(rho = ifelse(cop =="indep", 0, rho))

ci_char_names <- 
  c("est", "d1", "d2")
tibble_cisTruth_spearman <- 
  as_tibble(as.data.frame(as.table(array_cis_full))) %>%
  rename(N_index = Var1, cop_index = Var2, rho_index = Var3, 
         beta_index = Var4, m_index = Var5, 
         method_index = Var6, 
         ci_index = Var7, value = Freq) %>%
  mutate(
    N = (seq(1,N))[N_index],
    cop = copVec[cop_index],
    rho = rhoVec[rho_index],
    beta = betaVec[beta_index], 
    m = mVec[m_index],
    method = method_names[method_index],
    ci_char = ci_char_names[ci_index]
  ) %>%
  select(N, cop, rho, beta, m, method, ci_char, value) %>% 
  filter(!(cop == "indep" & rho == 0.4)) %>% 
  mutate(rho = ifelse(cop == "indep", 0, rho)) %>% 
  pivot_wider(names_from = "ci_char") %>% 
  full_join(tib_truth)

tibble_ci_covs_spearman <- 
  tibble_cisTruth_spearman %>% 
  mutate(cov = ifelse(est + d1 <= truth & 2*
                        est + d2 >= truth, 1, 0)) %>% 
  group_by(cop, rho, beta, m, method) %>% 
  summarise(emp_cov = mean(cov))

tibble_ci_covs_spearman$rho <-
  as.factor(tibble_ci_covs_spearman$rho)
tibble_ci_covs_spearman$cop <- 
  factor(tibble_ci_covs_spearman$cop, 
         levels = 
           c("indep", "gauss", "tcop", "gumbel"))
tibble_ci_covs_spearman <- bind_cols(tibble_ci_covs_spearman,
                                     tibble(char = "emp_cov")) %>% 
  rename(value = emp_cov)
plotCiEmpCov_spearman <- 
  tibble_ci_covs_spearman %>% 
  filter(beta == 0.5, method != "sb", rho != 0.2) %>% 
  ggplot(
    aes(x = m, y = value, col = method, linetype = rho))+
  facet_grid(cols = vars(cop), labeller =
               labeller(cop = facet_labels_cop))+
  geom_line(linewidth = 1)+
  #scale_y_continuous(limits = c(0.865, 0.98))+
  geom_abline(intercept = 0.95, slope = 0, col = "black", linetype = "longdash")+
  scale_color_manual(values = ownPalette)+
  scale_linetype_manual(values = c("dashed", "solid", "dotdash")) +
  scale_y_continuous(limits = c(0.75, 1))+
  themePlot+
  labs(#title = "Bootstrap confidence intervals",
    y = "Empirical Coverage", x = "Effective sample size m",
    col = "Confidence Interval:", linetype = expression(rho))
plotCiEmpCov_spearman 


#plot both relative avg width and empirical coverage
tibble_rel_avg_ci_emp_cov_spearman <- 
  bind_rows(tibble_rel_avg_ci_width_spearman, tibble_ci_covs_spearman)

facet_labels_chars <- as_labeller(
  c(
    "emp_cov" = "Empirical Coverage", 
    "rel_avg_width" = "Relative CI width"
  )
)


plotRelCiWidthCov_spearman <-
  tibble_rel_avg_ci_emp_cov_spearman %>% 
  filter(beta == 0, method != "sb", rho != 0.2) %>% 
  select(-rho) %>% 
  ggplot(aes(x = m, y = value, col = method)) +
  facet_grid(cols = vars(cop), rows = vars(char), 
             labeller = labeller(cop = facet_labels_cop, char = facet_labels_chars)) +
  
  # geom_abline für "emp_cov" nur in der ersten Zeile
  geom_abline(data = tibble_rel_avg_ci_emp_cov_spearman %>% 
                filter(char == "emp_cov") %>% mutate(dummy_x = NA, dummy_y = NA), 
              aes(intercept = 0.95, slope = 0), color = "black", linetype = "longdash") +
  
  # geom_abline für "rel_avg_width" nur in der zweiten Zeile
  geom_abline(data = tibble_rel_avg_ci_emp_cov_spearman %>% 
                filter(char == "rel_avg_width") %>% mutate(dummy_x = NA, dummy_y = NA), 
              aes(intercept = 1, slope = 0), color = color_db, linewidth = 1.1) +
  
  geom_line(linewidth = 1) +
  scale_color_manual(values = ownPalette) +
  scale_linetype_manual(values = c("dashed", "solid", "dotdash")) +
  themePlot +
  labs(y = "", x = "Effective sample size m", col = "Confidence Interval:", linetype = expression(rho))+
  theme(legend.position = "bottom")
plotRelCiWidthCov_spearman
if(TRUE){
  ggsave(plotRelCiWidthCov_spearman, 
         filename = "Figure F.7.pdf",
         device = "pdf", path = plotPath,
         width = 13, height = 7)
}
