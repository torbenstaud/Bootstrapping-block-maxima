library(ggplot2)
library(tidyverse)
library(latex2exp)
library(here)
library(ggpubr)
load(here("scripts/mean/bst/data/bstMeanTib"))
load(here("scripts/mean/bst/generateTrueVars/data/trueVarsTib"))
#load the tibbles and choose models----
bstMeanVarTib <- bstMeanTib
bstPlotTib <- bstMeanVarTib
bstPlotTib <- bstPlotTib %>% 
  left_join(trueVarsTib, by = c("gamma", "tsModel", "marginal", "m", 
                                       "r", "beta")) %>% 
  filter(
  gamma %in% c(-0.2, 0, 0.2), marginal == 1, k != 1, beta %in% c(0,0.5)
)
#create factors and plot----
##themes----
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

labsPlot <- labs(
  x = "Effective sample size m",
  col = "Bootstrap")
##factors----
bstPlotTib$k <- factor(bstPlotTib$k, levels = c(1,2,3,0),
                           labels = c(parse(text = TeX("db")), 
                                      parse(text = TeX("cb(2)")), 
                                      parse(text = TeX("cb(3)")), 
                                      parse(text = TeX("sb"))
                           )
)
bstPlotTib$marginal <- factor(bstPlotTib$marginal,
                                  levels = c(1,5),
                                  labels = c("GPD", "GEV")
)

bstPlotTib$beta <- factor(bstPlotTib$beta, 
                              levels = c(0, 0.25, 0.5, 0.75),
                              labels = c(parse(text = TeX("i.i.d.")), 
                                         parse(text = TeX("ARMAX $0.25$")), 
                                         parse(text = TeX("ARMAX $0.5$")), 
                                         parse(text = TeX("ARMAX $0.75$"))
                              )
)

bstPlotTib$gamma <- factor(bstPlotTib$gamma, 
                               levels = c(-0.2, -0.1, 0, 0.1, 0.2), 
                               labels = c(parse(text = TeX("$\\gamma = -0.2$")),
                                          parse(text = TeX("$\\gamma = -0.1$")),
                                          parse(text = TeX("$\\gamma = 0$")),
                                          parse(text = TeX("$\\gamma = 0.1$")),
                                          parse(text = TeX("$\\gamma = 0.2$"))
                                          ))



rescale <- 100
plotMeanFixRVarEst <- 
bstPlotTib %>% ggplot(aes(x = m, y = varBst*rescale, col = k))+
  geom_line()+ 
  geom_line(aes(x = m, y = trueVar * rescale), 
            col = "black", linetype = "dashed")+
  facet_grid(gamma~marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Variance * 100"),
        title = "Bootstrapping the variance of the mean estimator (r = 90 fixed)",
        linetype = "bla")+
  themePlot+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plotMeanFixRVarEst
if(FALSE){
  ggsave(plotMeanFixRVarEst, path = here("results/"), 
         filename = "plotMeanFixRVarEstMain.pdf", device = "pdf",
         width = 10, height = 10
  )
}
##relative plot----
plotRelTib <- bstPlotTib %>% mutate(relVar = varBst / trueVar)

plotMeanFixRVarEstRel <- 
plotRelTib %>% ggplot(aes(x = m, y = relVar, col = k))+
  geom_line()+ 
  geom_hline(yintercept = 1, col = "black", linetype = "dashed")+
  facet_grid(gamma~marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  scale_y_continuous(limits = c(0.6, 1.1))+
  labs( y = paste0("Relative variance"),
        title = "Bootstrapping the variance of the mean estimator (r = 90 fixed)",
        linetype = "bla")+
  themePlot+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plotMeanFixRVarEstRel
if(FALSE){
  ggsave(plotMeanFixRVarEstRel, path = here("results/"), 
         filename = "plotMeanBstVarRel.pdf", device = "pdf",
         width = 10, height = 10
  )
}

#Confidence intervals----
load(here("scripts/mean/bst/data/ciMeanTibRfix"))
ciMeanTibRfix

#create tibbles----

#create factors----
ciMeanTibRfix <- ciMeanTibRfix %>% filter(
  beta %in% c(0, 0.5)
)

ciMeanTibRfix <- ciMeanTibRfix %>% filter(
  gamma %in% c(-0.2,0, 0.2)
)


ciMeanTibRfix$k <- factor(ciMeanTibRfix$k, levels = c(2,3,0,1),
                             labels = c(parse(text = TeX("cb(2)")), 
                                        parse(text = TeX("cb(3)")), 
                                        parse(text = TeX("sb")),
                                        parse(text = TeX("db"))
                             )
)
ciMeanTibRfix$marginal <- factor(ciMeanTibRfix$marginal,
                                    levels = c(1,5),
                                    labels = c("GPD", "GEV")
)

ciMeanTibRfix$beta <- factor(ciMeanTibRfix$beta, 
                                levels = c(0, 0.25, 0.5, 0.75),
                                labels = c(parse(text = TeX("i.i.d.")), 
                                           parse(text = TeX("ARMAX $0.25$")), 
                                           parse(text = TeX("ARMAX $0.5$")), 
                                           parse(text = TeX("ARMAX $0.75$"))
                                )
)

ciMeanTibRfix$gamma <- factor(ciMeanTibRfix$gamma, 
                                 levels = c(0-0.2, -0.1, 0, 0.1, 0.2), 
                                 labels = c(parse(text = TeX("$\\gamma = -0.2$")),
                                            parse(text = TeX("$\\gamma = -0.1$")),
                                            parse(text = TeX("$\\gamma = 0$")),
                                            parse(text = TeX("$\\gamma = 0.1$")),
                                            parse(text = TeX("$\\gamma = 0.2$"))
                                            )
)
##plots----
ciCovPlot <- 
  ciMeanTibRfix %>%   
  filter(k != "sb", marginal !="GEV") %>% 
  ggplot()+
  geom_line(aes(x = m, y = empCov, col = k))+
  geom_hline(yintercept = 0.95, col = "black", linetype = "dashed")+
  facet_grid(gamma~marginal+beta,  labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Average empirical coverage"),
        title = "Confidence intervals for the mean")+
  scale_y_continuous(
    breaks = c(0.91, 0.93, 0.95)
  )+
  themePlot+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
ciCovPlot

#rescale <- 10
ciWidthPlot <-
  ciMeanTibRfix %>% 
  filter(k != "sb", marginal != "GEV") %>% 
  ggplot()+
  geom_line(aes(x = m, y = avgWidth, col = k))+
  facet_grid(gamma~beta+marginal, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Width"),
        title = "Average width of confidence intervals")+
  themePlot
ciWidthPlot

##now relative width Plot
ciWidthDb <- ciMeanTibRfix %>% filter(k == "db") %>% select(-c(k,empCov)) %>% 
  rename(dbWidth = "avgWidth")
ciWidthRel <- left_join(ciMeanTibRfix, ciWidthDb) %>% 
  filter(k %in% c("cb(2)", "cb(3)")) %>% mutate(relWidth = dbWidth/avgWidth)
ciWidthRelPlot <- ciWidthRel %>% filter(marginal != "GEV") %>% 
  ggplot()+
  geom_line(aes(x = m, y = relWidth, col = k))+
  geom_hline(yintercept = 1, col = "black", linetype = "dashed")+
  facet_grid(gamma~marginal+beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Relative average width"),
        title = "Average relative width of confidence intervals",
  )+
  scale_y_continuous(limits = c(1, 1.1), breaks = c(1,1.05, 1.1))+
  themePlot+
  theme(strip.text.x = element_blank(), plot.title = element_blank())

ciWidthRelPlot

combCiPlot <- 
  ggarrange(ciCovPlot, ciWidthRelPlot, nrow = 2, align = "v", common.legend = T, 
            legend = "bottom")
combCiPlot

if(FALSE){
  ggsave(combCiPlot, path = here("results/"), 
         filename = "plotMeanFixRCiMain.pdf", device = "pdf",
         width = 10, height = 10
  )
}