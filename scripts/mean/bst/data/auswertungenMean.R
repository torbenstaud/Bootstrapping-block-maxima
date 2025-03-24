library(here)
library(Rcpp)
library(dplyr)
library(tidyr)
library(ggplot2)
require(latex2exp)#for writing latex tau in the facets
library(ggpubr)
##plotting cosmetics
scaleFactor <- 1
textSize <- 20
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



#create QQ-Plots----
## source files: MACOS


## load data: name is "estArrayMean"----
#load("Z:/bootstrap/mean/data/r90/meanEst_gamma-0.2_Ts1_Distr1R.data")
#Parameters: gamma = -0.2, marginal = gpd, tsmod = armax, beta = 0.5, r = 90, m = 80

load(here("scripts/mean/bst/data/MeanEst_R90gamma-0.2_Ts3_Distr1_beta0.5.Rdata")) #meanArray
load(here("scripts/mean/bst/data/MeanBst_R90gamma-0.2_Ts3_Distr1_beta0.5_M80.Rdata"))  #estarray
sourceCpp(here("src/0sourceC.cpp"))

source(here("scripts/mean/bst/0parametersBstRFix.R"))

#N = 5*10^3
#kVec = c(0, 1,2,3)
#B = 10^3
#mVec = seq(40,100, by = 10)

#QQ-Plots and Histogramms----
#dims = N x mVec x kVec x B


# parameters: 90 (blocksize), 1(gpd), 3(tsmod), -0.2(gamma), 10^6(repetitions),
# 0.5(beta)

estsArray <- meanArray[,5,] #5 is the index for m = 80
bstArray <- estArray[,,]
truthMean <- getTruthC(90, 1,3,-0.2,10^6,0.5)
estsErrArray <- estsArray - truthMean
##create bootstrap error array
tmpJoinArray <- array(dim = c(N, length(kVec), B+1))
tmpJoinArray[,1,2] <- estsArray[,2]
tmpJoinArray[,c(1,3,4),1] <- estsArray[,1]
tmpJoinArray[,,seq(2,B+1)] <- bstArray
errorFun <- function(xx){
  xx[-1] - xx[1]
}
bstsErrArray <- tmpJoinArray %>% apply(c(1,2), errorFun)
#plot from the introduction
method_names <- c("sb", "db", "cb(2)", "cb(3)")
tib_estErr <- 
as_tibble(as.data.frame(as.table(estsErrArray))) %>%
  rename(N_index = Var1, 
         method_index = Var2, 
         value = Freq) %>% 
  mutate(N = N_index, method = method_names[method_index], 
         char = "est_err", gamma = -0.2, parameter = "mean") %>% 
  select(N, gamma, parameter, method, char, value)

tib_bstErr <- 
  as_tibble(as.data.frame(as.table(matrix(bstsErrArray, ncol = 4)))) %>%
  rename(N_index = Var1, 
         method_index = Var2, 
         value = Freq) %>% 
  mutate(N = N_index, method = method_names[method_index], 
         char = "bst_err", gamma = -0.2, parameter = "mean") %>% 
  select(N, gamma, parameter, method, char, value)

tib_fullErr <- bind_rows(tib_estErr, tib_bstErr)


facet_labels_method <- as_labeller(
  c(
    "db" = "db",
    "sb" = "Naive Sliding",
    "cb(2)" = "Circular"
  )
)
tib_fullErr$method <- factor(tib_fullErr$method, 
                             levels = c("sb", "cb(2)", "cb(3)", "db"))
tib_sbvscb <- 
  bind_rows(tib_fullErr %>% filter(method == "sb", char == "est_err"),
            tib_fullErr %>% filter(method == "cb(2)", char == "bst_err"))
plotIntroInc <- 
tib_fullErr %>% 
  filter(method != "db", method != "cb(3)") %>% 
  ggplot(aes(x = value, y= after_stat(density)))+
  geom_histogram(aes(x = value, fill = char), 
                 alpha = 0.4, position = "identity", bins = 30)+
  geom_density(aes(x = value, linetype = char), linewidth = 1)+
  facet_grid(cols = vars(method), labeller = labeller(method = facet_labels_method))+
  xlim(-0.2, 0.2)+
  scale_y_continuous(limits = c(0,11), breaks = c(0, 3.0, 6.0, 9.0))+
  labs(
    y = "Density",
    fill = "Error Type", linetype = "Error Type"
  )+
  scale_fill_manual(labels = c("Bootstrap", "Estimation"), values = c("#F8766D", "#00BFC4"))+
  scale_linetype_manual(labels = c("Bootstrap", "Estimation"), values = c("solid", "dashed"))+
  themePlot+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()
  )
plotIntroInc
if(F){
  
  ggsave(plotIntroInc, filename = sprintf("plotIncGpdMean-0.1Armax(0.5)Hists.pdf"),
         device = "pdf", path = here("results/"),
         width = 12.5, height = 4)
}


#for the supplement
##create bootstrap error array
tmpJoinArray <- array(dim = c(N, length(kVec), B+1))
tmpJoinArray[,,1] <- estsArray
tmpJoinArray[,,seq(2,B+1)] <- bstArray
errorFun <- function(xx){
  xx[-1] - xx[1]
}
bstsErrArray <- tmpJoinArray %>% apply(c(1,2), errorFun)



##function for qq-plot
createQQMean <- function(kEst, kBst, noY = F){
  #kEst and kBst in index form: 1 = sliding, 2 = disjoint, 3 = 2-max, 4 = 3-max 
  estErrors <- 10*quantile(estsErrArray[, kEst], 
                           seq(0.001,0.999, length.out = N))
  bstErrors <- 10*quantile(c(bstsErrArray[,,kBst]), 
                           seq(0.001,0.999, length.out = N))
  ggplot(mapping = aes(x = estErrors, y = bstErrors))+
    geom_point()+
    xlim(-1.5, 1.5)+
    scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(0, 1))+
    geom_abline(col = "red")+
    labs(
      y = "Quantiles * 10", 
      title = "QQ-Plot of the estimation error vs bootstrap error for the mean"
    )+
    themePlot+
    theme(plot.title = element_blank(),
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
createHistMean <- function(kEst, kBst, noY = F){
  kNames <- 
    c("sb", "db", "cb", "cb(3)")
  histTib <- tibble(type = rep(
    c(paste0("Estimation Error"), 
      paste0("Bootstrap Error")),
    c(N, N*B)),
    val = c(estsErrArray[,kEst], c(bstsErrArray[,,kBst])))
  histTib %>% ggplot(aes(x = val))+
    geom_histogram(aes(x = val, fill = type, y = after_stat(density)), alpha = 0.3, 
                   position = "identity")+
    geom_density(aes(x = val, linetype = type))+
    xlim(-0.2,0.2)+
    scale_y_continuous(limits = c(0,11), breaks = c(0, 3.0, 6.0, 9.0))+
    labs(
      x = paste0(kNames[kEst], " vs ", kNames[kBst]),
      y = "Density",
      fill = "Error Type:", linetype = "Error Type:"
    )+
    scale_linetype_manual(values = c("solid", "dashed"))+
    theme(plot.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
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
pH2 <- createHistMean(1,3, noY = T)
pH3 <- createHistMean(1,4, noY = T)
pH4 <- createHistMean(2,2, noY = T)

#create qq plots
pQ1 <- createQQMean(1,1)
pQ2 <- createQQMean(1,3, noY = T)
#pQ3 <- createQQMean(1,4, noY = T)
pQ4 <- createQQMean(2,2, noY = T)
#combine both qq and histogram plots
histQQPlots <- ggarrange(pQ1, pQ2, pQ4, pH1,pH2,pH4, 
                       ncol = 3, nrow = 2, common.legend = T, 
                       legend = "right", widths = c(1.3,1,1,1,1,1),
                       heights = c(1,1.35))
#pubhistQQPlots <- annotate_figure(histQQPlots, 
#                top = text_grob("Mean estimation error", 
#                                size = textSize, face = "bold"))

#pubhistQQPlots
histQQPlots
if(F){
  
ggsave(histQQPlots, filename = "plotHistQQMean.pdf",
       device = "pdf", path = here("results/"),
       width = 13, height = 7)
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

