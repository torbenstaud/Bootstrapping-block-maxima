library(here)
library(ggplot2)
library(tidyverse)
require(latex2exp)#for writing latex tau in the facets
library(ggpubr)
#-----------------------------------------------------------------------------
load(here("scripts/mle/fixedN/bst/generateTrueVars/data/trueVarsTib"))
load(here("scripts/mle/fixedN/bst/data/bstFreVarTib"))
source(here("src/theme.R"))
#create tibbles----

bstFreVarTibPl <- bstFreVarTib %>% 
  left_join(trueVarsTib) %>% 
  filter( beta %in% c(0, 0.5), k != 1) %>% 
  bind_cols(tibble(char = "sb"))


## for disjoint blocks----
load(here("scripts/mle/fixedN/bst/generateTrueVars/data/trueVarsTibDb"))#db true vars
bstFreVarTibPlDb <- bstFreVarTib %>% filter(k == 1) %>% 
  left_join(trueVarsTib, by = c("alpha", "tsModel", "marginal", "m", 
                                "r", "beta")) %>% 
  filter(
    alpha %in% c(0.5, 1, 1.5), beta %in% c(0,0.5)
  )%>% bind_cols(tibble(char = "db"))
bstFreVarTibPl <- bind_rows(bstFreVarTibPl, bstFreVarTibPlDb)
#create factors and plots----
##factors----
bstFreVarTibPl$k <- factor(bstFreVarTibPl$k, levels = c(0,2,3,1),
                           labels = c(parse(text = TeX("sb")), 
                                      parse(text = TeX("cb(2)")), 
                                      parse(text = TeX("cb(3)")), 
                                      parse(text = TeX("db"))
                           )
)
ownPalette <- #based on dark2
  c("cb(2)" = "#F8766D",  
    "cb(3)" = "#7CAE00",  
    "db" = "#00BFC4",  
    "sb" = "#C77CFF")

bstFreVarTibPl$marginal <- factor(bstFreVarTibPl$marginal,
                                  levels = c(3,4),
                                  labels = c("Fréchet", "Pareto")
)

bstFreVarTibPl$beta <- factor(bstFreVarTibPl$beta, 
                              levels = c(0, 0.25, 0.5, 0.75),
                              labels = c(parse(text = TeX("i.i.d.")), 
                                         parse(text = TeX("ARMAX $0.25$")), 
                                         parse(text = TeX("ARMAX $0.5$")), 
                                         parse(text = TeX("ARMAX $0.75$"))
                              )
)

bstFreVarTibPl$alpha <- factor(bstFreVarTibPl$alpha, 
                               levels = c(0.5, 1, 1.5), 
                               labels = c(parse(text = TeX("$\\alpha_0 = 0.5$")),
                                          parse(text = TeX("$\\alpha_0 = 1$")),
                                          parse(text = TeX("$\\alpha_0 = 1.5$")))
)
bstFreVarTibPl$char <- factor(bstFreVarTibPl$char, 
                              levels = c("sb", "db"), 
                              labels = c(parse(text = TeX("$\\sigma^2_{sb}$")),
                                         parse(text = TeX("$\\sigma^2_{db}$"))
                              ))
##plots----
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
labsPlot <- labs(
  x = "Effective sample size m",
  col = "Bootstrap:")
rescale <- 100

bstVarPlot <- 
  bstFreVarTibPl %>% mutate(asVar = m*varBst, asTrueVar = m*trueVar) %>% 
  filter(alpha != 2) %>% 
  ggplot()+
  geom_line(aes(x = m, y = varBst*rescale, col = k),
            linewidth = 1.1)+
  geom_line(aes(x = m, y = trueVar*rescale), col = "black", linetype = "dashed",
            linewidth = 1.1)+
  facet_grid(alpha~char+marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Variance * ",rescale),
        title = "Bootstrapping the variance of the shape estimator")+
  themePlot+
  scale_color_manual(values = ownPalette)+
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    plot.title = element_blank()
  )+
  scale_x_continuous(
    breaks = c(25, 50, 75)
  )
bstVarPlot
if(FALSE){
  ggsave(bstVarPlot, path = here("results/"), 
         filename = "plotFreFixNVarEstSupp.pdf", device = "pdf",
         width = 15, height = 10, 
  )
}
#now for the main paper: no Fréchet marginals-----

bstVarPlotMain <- 
  bstFreVarTibPl %>% filter(marginal == "Pareto") %>% 
  mutate(asVar = m*varBst, asTrueVar = m*trueVar) %>% 
  filter(alpha != 2) %>% 
  ggplot()+
  geom_line(aes(x = m, y = varBst*rescale, col = k),
            linewidth = 1.1)+
  geom_line(aes(x = m, y = trueVar*rescale), col = "black", linetype = "dashed",
            linewidth = 1.1)+
  facet_grid(alpha~char+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Variance * ",rescale),
        title = "Bootstrapping the variance of the shape estimator")+
  themePlot+
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "right",
    plot.title = element_blank()
  )+
  scale_color_manual(values = ownPalette)+
  scale_x_continuous(
    breaks = c(25, 50, 75)
  )
bstVarPlotMain
if(FALSE){
  ggsave(bstVarPlotMain, path = here("results/"), 
         filename = "plotFreFixNVarEstMain.pdf", device = "pdf",
         width = 10, height = 10, 
  )
}
#for jrssb
if(FALSE){
  ggsave(bstVarPlotMain, path = here("results/"), 
         filename = "plotFreFixNVarEstMainJrssb.pdf", device = "pdf",
         width = 15, height = 10, 
  )
}


#Confidence Intervals----
load(here("scripts/mle/fixedN/bst/data/CiFreNfix"))
ciFreTibNfix
##sizes----
textSize <- 15
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
##create factors----
ciFreTibNFixPlot <- ciFreTibNfix %>% filter(
  beta %in% c(0, 0.5), 
  m >=25, m<=83
)




ciFreTibNFixPlot$k <- factor(ciFreTibNFixPlot$k, levels = c(0,2,3,1),
                             labels = c(parse(text = TeX("sb")),
                                        parse(text = TeX("cb(2)")), 
                                        parse(text = TeX("cb(3)")), 
                                        parse(text = TeX("db"))
                             )
)
ciFreTibNFixPlot$marginal <- factor(ciFreTibNFixPlot$marginal,
                                    levels = c(3,4),
                                    labels = c("Fréchet", "Pareto")
)

ciFreTibNFixPlot$beta <- factor(ciFreTibNFixPlot$beta, 
                                levels = c(0, 0.25, 0.5, 0.75),
                                labels = c(parse(text = TeX("i.i.d.")), 
                                           parse(text = TeX("ARMAX $0.25$")), 
                                           parse(text = TeX("ARMAX $0.5$")), 
                                           parse(text = TeX("ARMAX $0.75$"))
                                )
)

ciFreTibNFixPlot$alpha <- factor(ciFreTibNFixPlot$alpha, 
                                 levels = c(0.5, 1, 1.5), 
                                 labels = c(parse(text = TeX("$\\alpha_0 = 0.5$")),
                                            parse(text = TeX("$\\alpha_0 = 1$")),
                                            parse(text = TeX("$\\alpha_0 = 1.5$")))
)
##plots----
ciCovPlot <- 
  ciFreTibNFixPlot %>%  filter( m != 31) %>% 
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = empCov, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 0.95, col = "black", linetype = "dashed",
             linewidth = 1.1)+
  facet_grid(alpha~marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Empirical coverage"),
        title = "Confidence intervals for the shape")+
  themePlot+
  scale_color_manual(values = ownPalette)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "right"
  )+
  scale_y_continuous(
    breaks = c(0.85, 0.9, 0.95)
  )
ciCovPlot

#rescale <- 10
ciWidthPlot <-
  ciFreTibNFixPlot %>%  
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = avgWidth, col = k),
            linewidth = 1.1)+
  facet_grid(alpha~marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Width"),
        title = "Average width of confidence intervals")+
  scale_color_manual(values = ownPalette)+
  themePlot+
  theme(
    legend.position = "right"
  )
ciWidthPlot

##now relative width Plot
ciWidthDb <- ciFreTibNFixPlot %>% filter(k == "db") %>% select(-c(k,empCov)) %>% 
  rename(dbWidth = "avgWidth")
ciWidthRel <- left_join(ciFreTibNFixPlot, ciWidthDb) %>% 
  filter(k %in% c("cb(2)", "cb(3)")) %>% mutate(relWidth = dbWidth/avgWidth)
ciWidthRelPlot <- ciWidthRel %>% filter( m != 31) %>% 
  ggplot()+
  geom_line(aes(x = m, y = relWidth, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 1, col = "#00BFC4", linetype = "dashed",
             linewidth = 1.1)+
  facet_grid(alpha~marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Relative average width"),
        title = "Average relative width of confidence intervals",
  )+
  scale_y_continuous(limits = c(0.95, 1.25), breaks = c(1,1.1,1.2))+
  scale_x_continuous(breaks = c(25, 50, 75))+
  themePlot+
  scale_color_manual(values = ownPalette)+
  theme(strip.text.x = element_blank(), plot.title = element_blank(), 
        axis.text.x = element_text(angle = 90))

ciWidthRelPlot

combCiPlot <- 
  ggarrange(ciCovPlot, ciWidthRelPlot, nrow = 2, align = "v", common.legend = T, 
            legend = "right",
            heights = c(1, 1.05))
combCiPlot

if(FALSE){
  ggsave(combCiPlot, path = here("results/"), 
         filename = "plotFreFixNCiSupp.pdf", device = "pdf",
         width = 12, height = 12, 
  )
}

#now for the main paper, where we omit the Fréchet marginal case
ciFreTibNFixPlotMain <- ciFreTibNFixPlot %>% filter( m != 31, marginal == "Pareto")
ciCovPlotMain <- 
  ciFreTibNFixPlotMain %>%  filter( m != 31, marginal == "Pareto") %>% 
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = empCov, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 0.95, col = "black", linetype = "dashed",
             linewidth = 1.1)+
  facet_grid(alpha~beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Empirical coverage"),
        title = "Confidence intervals for the shape")+
  themePlot+
  scale_color_manual(values = ownPalette)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank())+
  scale_y_continuous(
    breaks = c(0.85, 0.9, 0.95)
  )
ciCovPlotMain

#rescale <- 10
ciWidthPlotMain <-
  ciFreTibNFixPlotMain %>%  filter( m != 31, marginal == "Pareto") %>%
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = avgWidth, col = k),
            linewidth = 1.1)+
  facet_grid(alpha~marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Width"),
        title = "Average width of confidence intervals")+
  themePlot+
  scale_color_manual(values = ownPalette)
ciWidthPlotMain

##now relative width Plot
ciWidthDbMain <- ciFreTibNFixPlotMain %>% filter(k == "db") %>% select(-c(k,empCov)) %>% 
  rename(dbWidth = "avgWidth")
ciWidthRelMain <- left_join(ciFreTibNFixPlotMain, ciWidthDbMain) %>% 
  filter(k %in% c("cb(2)", "cb(3)")) %>% mutate(relWidth = dbWidth/avgWidth)
ciWidthRelPlotMain <- ciWidthRelMain %>% filter( m != 31) %>% 
  ggplot()+
  geom_line(aes(x = m, y = relWidth, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 1, col = "#00BFC4", linetype = "dashed",
             linewidth = 1.1)+
  facet_grid(alpha~beta+marginal, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Relative average width"),
        title = "Average relative width of confidence intervals",
  )+
  scale_y_continuous(limits = c(0.95, 1.25), breaks = c(1,1.1,1.2))+
  scale_x_continuous(breaks = c(25, 50, 75))+
  themePlot+
  scale_color_manual(values = ownPalette)+
  theme(strip.text.x = element_blank(), plot.title = element_blank(), 
        axis.text.x = element_text(angle = 90))

ciWidthRelPlotMain

combCiPlotMain <- 
  ggarrange(ciCovPlotMain, ciWidthRelPlotMain, nrow = 2, align = "v", common.legend = T, 
            legend = "right",
            heights = c(1, 1.05))
combCiPlotMain

if(FALSE){
  ggsave(combCiPlotMain, path = here("results/"), 
         filename = "plotFreFixNCiMain.pdf", device = "pdf",
         width = 10, height = 10, 
  )
}
#for jrssb
ciCovPlotMainJrssb <- 
  ciFreTibNFixPlotMain %>%  filter( m != 31, marginal == "Pareto") %>% 
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = empCov, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 0.95, col = "black", linetype = "dashed",
             linewidth = 1.1)+
  facet_grid(alpha~beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Empirical coverage"),
        title = "Confidence intervals for the shape")+
  themePlot+
  scale_color_manual(values = ownPalette)+
  theme(strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_blank()
        )+
  scale_y_continuous(
    breaks = c(0.85, 0.9, 0.95)
  )
ciCovPlotMainJrssb


ciWidthRelPlotMainJrssb <- ciWidthRelMain %>% filter( m != 31) %>% 
  ggplot()+
  geom_line(aes(x = m, y = relWidth, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 1, col = "#00BFC4", linetype = "dashed",
             linewidth = 1.1)+
  facet_grid(alpha~beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Relative average width"),
        title = "Average relative width of confidence intervals",
  )+
  scale_y_continuous(limits = c(0.95, 1.25), breaks = c(1,1.1,1.2))+
  scale_x_continuous(breaks = c(25, 50, 75))+
  themePlot+
  scale_color_manual(values = ownPalette)+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ciWidthRelPlotMainJrssb

combCiPlotMainJrssb <- 
  ggarrange(ciCovPlotMainJrssb, ciWidthRelPlotMainJrssb, ncol = 2, align = "h", 
            common.legend = T, 
            legend = "bottom")
combCiPlotMainJrssb

if(FALSE){
  ggsave(combCiPlotMainJrssb, path = here("results/"), 
         filename = "plotFreFixNCiMainJrssb.pdf", device = "pdf",
         width = 13, height = 7
  )
}