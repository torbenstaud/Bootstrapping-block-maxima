library(here)
library(tidyverse)
textSize <- 20
source(here("backend/0plot_theme.R"))
plotPath <- file.path("Supplement", "plots")

dataPath <- here(file.path("backend", "data", "Supplement", "F.2."))



#Alle harcoded variablen definieren
copVec <- c("indep", "gauss", "gumbel", "tcop")
betaVec <- c(0, 0.5)
rhoVec <- c(0.2, 0.4) #der indep Fall ist besonders, da kt hier 0 ist
mVec <- seq(40,100, by = 10)
r <- 365
N <- 10**3
B <- 10**3
NSeedChunkSize <- 500
NSeedChunkN <- N/NSeedChunkSize
NSeedChunkVec <- seq(1, NSeedChunkN)
Time <- 100
level <- 0.05

#array mit CIs laden und zum tibble trafon
load(here(dataPath, "tib_ests_full")) 
tib_ests_full <- tib_ests_full %>% 
  filter(!(rho == 0.4 & cop == "indep")) %>% 
  mutate(method = ifelse(method == "sb-cb", "cb", method))


#Wahrheiten als tibble laden
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
tib_truth <- tib_truthTemp%>% 
  rename(rho = rho_innov, cop = coptype) %>% filter(!(cop == "indep" & rho == 0.4))


tib_ests <- full_join(tib_ests_full, tib_truth)

tib_ests_chars_full <-
  tib_ests %>% mutate(bias = value - rho_bm) %>% 
  group_by(cop, rho, beta, m, method) %>% 
  summarise(emp_mean = mean(value), truth = mean(rho_bm), variance = var(value),
            bias = truth - emp_mean, mse = variance + bias**2, bias_sq = bias**2, 
            bias_par = bias**2/mse)
tib_relmse <- 
  tib_ests_chars_full %>% 
  select(c(cop,rho,beta,m,method,mse)) %>% rename(value = mse) %>% 
  pivot_wider(names_from = method) %>% 
  mutate(sb = sb/db, `cb` = `cb`/db, db = 1) %>% 
  pivot_longer(cols = c(db, sb, `cb`), names_to = "method", values_to = "rel_mse")

tib_ests_chars_full <- full_join(tib_ests_chars_full, tib_relmse)




########


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
tib_ests_chars_full$cop <- factor(tib_ests_chars_full$cop, 
                                  levels = 
                                    c("indep", "gauss", "tcop", "gumbel"))

facet_labels_name <-
  as_labeller(
    c("bias_par" = "Bias Contribution",
      "rel_mse" = "Relative MSE",
      "mse" = "MSE"))

plotEstsMse_spearman <- 
  tib_ests_chars_full %>% ungroup %>% 
  mutate(rho = ifelse(cop == "indep" & rho == 0.2, 0, rho)) %>% 
  filter(beta == 0, rho != 0.2) %>%
  select(c(cop, m, method, mse, rel_mse)) %>% 
  pivot_longer(cols = c(mse, rel_mse)) %>% 
  ggplot(aes(x=m, y = value, col = method)) +
  # Unterschiedliche Linewidths für bias_sq und variance im varbias-Facet
  geom_line(linewidth = 1.1) +
  
  # Anpassung der Farben und Linien
  scale_color_manual(values = ownPalette, na.translate = FALSE) +
  scale_linetype_manual(values = c("dashed", "solid", "dotdash"), 
                        na.translate = FALSE) +
  
  # Individuelle Achsenlimits pro Facet (NUR für bias_par)
  facet_grid(rows = vars(name), cols = vars(cop), scales = "free_y",
             labeller = labeller(cop = facet_labels_cop,
                                 name = facet_labels_name)) +
  
  # Achsenbeschriftungen und Legenden
  labs(
    x = "Effective sample size m",
    y = "Value",
    linetype = expression(rho),
    col = "Estimator:"
  ) +
  
  # Design-Anpassungen
  themePlot+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_blank()
  )
plotEstsMse_spearman
if(T){
  
  ggsave(plotEstsMse_spearman, 
         filename = "Figure F.6.pdf",
         device = "pdf", path = plotPath,
         width = 13, height = 7)
  
}
