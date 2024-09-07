library(here)
library(Rcpp)
library(dplyr)
library(tidyr)
library(ggplot2)
require(latex2exp)#for writing latex tau in the facets
library(ggpubr)
library(RColorBrewer)
load(here("scripts/mean/est/data/fullVarTibMeanEst")) #fullVarTibMeanEst
fullVarTibMeanEst

fullVarTibMeanEst <- fullVarTibMeanEst %>% 
  filter(marginal == 1, 
         gamma %in% c(-0.2, 0, 0.2),
         beta %in% c(0, 0.5))
#create factors for plots
fullVarTibMeanEst$beta <- fullVarTibMeanEst$beta %>% factor(
  levels = c(0, 0.5), 
  labels = c("i.i.d.", 
             parse(text = TeX("ARMAX $0.5$")))
)
fullVarTibMeanEst$marginal <- factor(fullVarTibMeanEst$marginal,
                                     levels = c(1,5),
                                     labels = c("GPD", "GEV")
)

fullVarTibMeanEst$k <- fullVarTibMeanEst$k %>% factor(
  levels = c(0,2,3,1), 
  labels = c("sb", "cb(2)", "cb(3)", "db")
)
ownPalette <- #based on dark2
  c("cb(2)" = "#F8766D",  
    "cb(3)" = "#7CAE00",  
    "db" = "#00BFC4",  
    "sb" = "#C77CFF")

fullVarTibMeanEst$gamma <- fullVarTibMeanEst$gamma %>% factor(
  levels = c(-0.2, -0.1, 0, 0.1, 0.2),
  labels = c(parse(text = TeX("$\\gamma = -0.2$")),
             parse(text = TeX("$\\gamma = -0.1$")),
             parse(text = TeX("$\\gamma = 0$")),
             parse(text = TeX("$\\gamma = 0.1$")),
             parse(text = TeX("$\\gamma = 0.2$")))
)

textSize <- 20
themePlot <- theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2),
                   strip.background = element_rect(color = "black", 
                                                   fill = "lightgrey", size = 0.2),
                   axis.title.x = element_text(size = textSize),
                   axis.title.y = element_text(size = textSize),
                   axis.text.y =element_text(size=textSize), 
                   axis.text.x =element_text(size=textSize),
                   strip.text.x = element_text(size = textSize),
                   strip.text.y = element_text(size = textSize),
                   plot.title = element_text(hjust = 0.5, size = textSize, 
                                             face = "bold"), 
                   #panel.background = element_rect(rgb(0.95, 0.95, 0.95, alpha = 1)),
                   legend.position = "right",
                   legend.title = element_text(size = textSize),
                   legend.text = element_text(size = textSize))

labsPlot <- labs(title = "Mean estimation (r = 90 fixed)",
                 x = "Effective sample size m",
                 col = "Estimator:")


plotAbsMse <- 
  fullVarTibMeanEst %>% mutate(biasSq = biasEst^2) %>% 
  pivot_longer(cols = c("biasEst", "varEst", "mse", "biasSq"),
               names_to = "char", values_to = "val") %>% 
  filter(char == "mse") %>% 
  mutate(val = val * 10^2) %>% 
  ggplot(aes(x = m, y = val, col = k, linetype = char))+
  geom_line(linewidth = 1.1)+
  facet_grid(gamma~beta, scales = "free_y", labeller=label_parsed)+
  labsPlot + 
  labs(y = "MSE * 100")+
  themePlot+
  scale_color_manual(values = ownPalette)+
  guides(linetype = "none")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank())

plotAbsMse
# fullVarTibMeanEst %>% mutate(biasSq = biasEst^2) %>% 
#   pivot_longer(cols = c("biasEst", "varEst", "mse", "biasSq", "ratio"),
#                names_to = "char", values_to = "val") %>% 
#   filter(char == "ratio", k != "sb") %>% 
#   ggplot(aes(x = m, y = val, col = k, linetype = char))+
#   geom_line()+
#   facet_grid(gamma~tsMod, scales = "free_y", labeller=label_parsed)+
#   labsPlot + 
#   labs(y = "MSE * 100")+
#   themePlot+
#   guides(linetype = "none")
# now relatively
varTibDb <- fullVarTibMeanEst %>% filter(k == "db") %>% 
  select(gamma, tsMod, marginal, beta, m, r, mse) %>% 
  rename(dbMSE = "mse")
fullVarTibMeanEstRel <- left_join(fullVarTibMeanEst, varTibDb) %>% 
  filter(k != "db") %>% 
  mutate(ratioMSE = dbMSE / mse)

plotRelMse <- 
  fullVarTibMeanEstRel %>% 
  filter(tsMod != 10) %>% 
  ggplot(aes(x = m, y = ratioMSE, col = k))+
  geom_line(linewidth = 1.1)+
  facet_grid(gamma~marginal+beta, labeller = label_parsed)+
  labs( y = "Relative MSE")+
  ylim(1, 1.15)+
  scale_y_continuous(breaks = c(1, 1.05, 1.1))+
  labsPlot+
  geom_hline(yintercept = 1, color = "#00BFC4", linetype = "longdash",
             linewidth = 1.1)+
  themePlot+
  theme(
    plot.title = element_blank(),
    strip.background.x = element_blank(), 
    strip.text.x = element_blank()
  )+
  scale_color_manual(values = ownPalette)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plotRelMse
plotAbsRelMse <- ggarrange(plotAbsMse, plotRelMse, align = "v", ncol = 1,
                           common.legend = T, legend = "right",
                           heights = c(1, 1.15))
plotAbsRelMse

if(FALSE){
  ggsave(plotAbsRelMse, filename = "plotMeanAbsRelMse.pdf",
         device = "pdf", path = here("results/"),
         width = 10, height = 12)
}

# bias
plotRelMse2 <- 
  fullVarTibMeanEstRel %>% 
  filter(tsMod != 10) %>% 
  ggplot(aes(x = m, y = ratioMSE, col = k))+
  geom_line(linewidth = 1.1)+
  facet_grid(gamma~marginal+beta, labeller = label_parsed)+
  labs( y = "Relative MSE")+
  ylim(1, 1.15)+
  scale_y_continuous(breaks = c(1, 1.05, 1.1))+
  labsPlot+
  geom_hline(yintercept = 1, color = "black", linetype = "longdash",
             linewidth = 1.1)+
  themePlot+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
plotRelBias <- 
  fullVarTibMeanEst %>% mutate(biasEst = biasEst^2/mse) %>% 
  filter(k != "db") %>% 
  pivot_longer(cols = c("biasEst", "varEst", "mse"),
               names_to = "char", values_to = "val") %>% 
  filter(char == "biasEst") %>% 
  mutate(val = val * 100) %>% 
  ggplot(aes(x = m, y = val, col = k, linetype = char))+
  geom_line(linewidth = 1.1)+
  facet_grid(gamma~beta, labeller = label_parsed)+
  labs(y = "Bias^2/MSE * 100")+
  labsPlot +
  themePlot+
  theme(strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        plot.title =  element_blank())+
  #scale_x_continuous(breaks = c(50,60, 7080 , 100))+
  ylim(0,0.8)+
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6))+
  guides(linetype = "none")

plotRelMseBias <- ggarrange(plotRelMse2, plotRelBias, align = "v", ncol = 1,
                            common.legend = T, legend = "right",
                            heights = c(1, 1.15))

plotRelMseBias


if(FALSE){
  ggsave(plotRelMseBias, filename = "plotMeanRelMseBias.pdf",
         device = "pdf", path = here("results/"),
         width = 10, height = 12)
}

