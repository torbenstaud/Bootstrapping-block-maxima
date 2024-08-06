library(here)
library(Rcpp)
library(dplyr)
library(tidyr)
library(ggplot2)
require(latex2exp)#for writing latex tau in the facets
library(ggpubr)
#create QQ-Plots----
## source files: MACOS


## load data: name is "estArrayMean"----
#load("Z:/bootstrap/mean/data/r90/meanEst_gamma-0.2_Ts1_Distr1R.data")
load(here("scripts/mean/data/meanEst_gamma-0.2_Ts1_Distr1R.data")) 
source(here("scripts/mean/0parameters.R"))

N <- 10^3
#QQ-Plots and Histogramms----
#dims = N x mVec x kVec x B

gammaInd <- 1 #we want gamma = -0.2
gamma <- gammaVec[gammaInd]
errorFun <- function(xx){
  xx[-1] - xx[1]
}
truthMean <- getTruthC(90, 1,1,-0.2,10^6,0.5)
#explore
scaleFactor <- 2
textSize <- 12
themePlot <- theme(panel.border = element_rect(color = "black", fill = NA),
                   strip.background = element_rect(color = "black", 
                                                   fill = "lightgrey"),
                   axis.title.x = element_text(size = scaleFactor*textSize, 
                                               face = "plain"),
                   axis.title.y = element_text(size = scaleFactor*textSize, 
                                               face = "plain"),
                   axis.text.y =element_text(size = scaleFactor*textSize), 
                   axis.text.x =element_text(size = scaleFactor*textSize),
                   strip.text.x = element_text(size = scaleFactor*textSize),
                   strip.text.y = element_text(size = scaleFactor*textSize),
                   plot.title = element_text(hjust = 0.5, size = scaleFactor*textSize), 
                   #panel.background = element_rect(rgb(0.95, 0.95, 0.95, alpha = 1)),
                   legend.position = "right",
                   legend.title = element_text(size = scaleFactor*textSize),
                   legend.text = element_text(size = scaleFactor*textSize))

createQQMean <- function(kEst, kBst, mInd = 4, noY = F){
  #kEst and kBst in index form: 1 = sliding, 2 = disjoint, 3 = 2-max, 4 = 3-max 
  estErrors <- 10*quantile(estArrayMean[, mInd, kEst,1] - truthMean, 
                        seq(0.001,0.999, length.out = N))
  bstErrors <- 10*quantile(c(apply(estArrayMean[,mInd, kBst,], c(1), errorFun)), 
                        seq(0.001,0.999, length.out = N))
  ggplot(mapping = aes(x = estErrors, y = bstErrors))+
    geom_point()+
    xlim(-1.5, 1.5)+
    ylim(-1.5, 1.5)+
    geom_abline(col = "red")+
    labs(
      y = "Quantiles * 10", 
      title = "QQ-Plot of the estimation error vs bootstrap error for the mean"
    )+
    themePlot+
    theme(plot.title = element_blank(),
          axis.title.y = element_text(size = scaleFactor*textSize), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()
          )+
    if(noY == T){
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank())
    }
}
createHistMean <- function(kEst, kBst, mInd = 4, noY = F){
  kNames <- 
    c("sb", "db", "cb(2)", "cb(3)")
  estErrors <- estArrayMean[, mInd, kEst,1] - truthMean
  bstErrors <- c(apply(estArrayMean[,mInd, kBst,], c(1), errorFun))
  histTib <- tibble(type = rep(
    c(paste0("Estimation Error"), 
      paste0("Bootstrap Error")),
                               c(N, N*B)),
                    val = c(estErrors, bstErrors))
  histTib %>% ggplot(aes(x = val))+
    geom_histogram(aes(x = val, fill = type, y = after_stat(density)), alpha = 0.3, 
                   position = "identity")+
    geom_density(aes(x = val, linetype = type))+
    xlim(-0.2,0.2)+
    ylim(0, 11)+
    labs(
      x = paste0(kNames[kEst], " vs ", kNames[kBst]),
      y = "Density",
      fill = "Error Type", linetype = "Error Type"
    )+
    scale_linetype_manual(values = c("solid", "dashed"))+
    theme(plot.title = element_blank()
          )+
    themePlot+
    if(noY){
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    }
}


#create histogram plots
pH1 <- createHistMean(1,1)
pH2 <- createHistMean(1,4, noY = T)
pH3 <- createHistMean(4,4, noY = T)
pH4 <- createHistMean(2,2, noY = T)

#create qq plots
pQ1 <- createQQMean(1,1)
pQ2 <- createQQMean(1,4, noY = T)
pQ3 <- createQQMean(4,4, noY = T)
pQ4 <- createQQMean(2,2, noY = T)
#combine both qq and histogram plots
histQQPlots <- ggarrange(pQ1, pQ2, pQ3, pQ4, pH1,pH2,pH3,pH4, 
                       ncol = 4, nrow = 2, common.legend = T, 
                       legend = "bottom", widths = c(1,1), heights = c(1,1.2))
pubhistQQPlots <- annotate_figure(histQQPlots, 
                top = text_grob("Mean estimation (r = 90 fixed)", 
                                size = textSize*scaleFactor, face = "bold"))

pubhistQQPlots

if(F){
  
ggsave(pubhistQQPlots, filename = "plotHistQQMean.pdf",
       device = "pdf", path = here("results/"),
       scale = scaleFactor-0.3)
}

#vergleiche schätzer----
load(here("scripts/mean/data/processedData.RData"))

varTibMeanEst <- varTibMeanEst %>% 
  filter(marginal == 1, 
         gamma %in% c(-0.2, 0, 0.2))
#create factors for plots
varTibMeanEst$tsModel <- varTibMeanEst$tsModel %>% factor(
  levels = c(1,3), 
  labels = c("i.i.d.", 
             parse(text = TeX("ARMAX $0.5$")))
)
varTibMeanEst$k <- varTibMeanEst$k %>% factor(
  levels = c(1,2,3,0), 
  labels = c("db", "cb(2)", "cb(3)", "sb")
)

varTibMeanEst$gamma <- varTibMeanEst$gamma %>% factor(
  levels = c(-0.2, -0.1, 0, 0.1, 0.2),
  labels = c(parse(text = TeX("$\\gamma = -0.2$")),
             parse(text = TeX("$\\gamma = -0.1$")),
             parse(text = TeX("$\\gamma = 0$")),
             parse(text = TeX("$\\gamma = 0.1$")),
             parse(text = TeX("$\\gamma = 0.2$")))
)

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

labsPlot <- labs(title = "Mean estimation (r = 90 fixed)",
                 x = "Effective sample size m",
                 col = "Estimator")


varTibMeanEst %>% mutate(biasSq = empBias^2) %>% 
  pivot_longer(cols = c("empBias", "empVar", "empMSE", "biasSq"),
                              names_to = "char", values_to = "val") %>% 
  filter(char == "empMSE") %>% 
  mutate(val = val * 10^2) %>% 
  ggplot(aes(x = m, y = val, col = k, linetype = char))+
  geom_line()+
  facet_grid(gamma~tsModel, scales = "free_y", labeller=label_parsed)+
  labsPlot + 
  labs(y = "MSE * 100")+
  themePlot+
  guides(linetype = "none")
varTibMeanEst %>% mutate(biasSq = empBias^2) %>% 
  pivot_longer(cols = c("empBias", "empVar", "empMSE", "biasSq", "ratio"),
               names_to = "char", values_to = "val") %>% 
  filter(char == "ratio", k != "sb") %>% 
  ggplot(aes(x = m, y = val, col = k, linetype = char))+
  geom_line()+
  facet_grid(gamma~tsModel, scales = "free_y", labeller=label_parsed)+
  labsPlot + 
  labs(y = "MSE * 100")+
  themePlot+
  guides(linetype = "none")
# now relatively
varTibDb <- varTibMeanEst %>% filter(k == "db") %>% 
  select(gamma, tsModel, marginal, m, r, empMSE) %>% 
  rename(dbMSE = "empMSE")
varTibMeanEstRel <- left_join(varTibMeanEst, varTibDb) %>% 
  filter(k != "db") %>% 
  mutate(ratioMSE = dbMSE / empMSE)

plotRelMse <- 
  varTibMeanEstRel %>% 
    filter(tsModel != 10) %>% 
    ggplot(aes(x = m, y = ratioMSE, col = k))+
    geom_line()+
    facet_grid(gamma~tsModel, labeller = label_parsed)+
    labs( y = "Relative MSE")+
    ylim(1, 1.15)+
    scale_y_continuous(breaks = c(1, 1.05, 1.1))+
    labsPlot+
    geom_hline(yintercept = 1, color = "black", linetype = "dashed")+
    themePlot+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())

# bias
plotRelBias <- 
  varTibMeanEst %>% mutate(empBias = empBias^2/empMSE) %>% 
  pivot_longer(cols = c("empBias", "empVar", "empMSE"),
               names_to = "char", values_to = "val") %>% 
  filter(char == "empBias") %>% 
  mutate(val = val * 100) %>% 
  ggplot(aes(x = m, y = val, col = k, linetype = char))+
  geom_line()+
  facet_grid(gamma~tsModel, labeller = label_parsed)+
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

plotRelMseBias <- ggarrange(plotRelMse, plotRelBias, align = "v", ncol = 1,
          common.legend = T, legend = "bottom")

plotRelMseBias
if(FALSE){
  ggsave(plotRelMseBias, filename = "plotMeanRelMseBias.pdf",
         device = "pdf", path = here("results/"),
         width = 10, height = 10)
}

