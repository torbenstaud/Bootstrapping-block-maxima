library(ggplot2)
library(tidyverse)
library(latex2exp)
library(here)
library(ggpubr)
load(here("scripts/mean/bst/data/bstMeanTib"))
load(here("scripts/mean/bst/generateTrueVars/data/trueVarsTib"))#sb true vars


#load the tibbles and choose models----
bstMeanVarTib <- bstMeanTib
bstPlotTib <- bstMeanVarTib
bstPlotTib <- bstPlotTib %>% 
  left_join(trueVarsTib, by = c("gamma", "tsModel", "marginal", "m", 
                                "r", "beta")) %>% 
  filter(
    gamma %in% c(-0.2, 0, 0.2), marginal == 1, k != 1, beta %in% c(0,0.5)
  ) %>% bind_cols(tibble(char = "sb"))
## for disjoint blocks----
load(here("scripts/mean/bst/generateTrueVars/data/trueVarsTibDb"))#db true vars
bstPlotTibDb <- bstMeanVarTib %>% filter(k == 1) %>% 
  left_join(trueVarsTib, by = c("gamma", "tsModel", "marginal", "m", 
                                "r", "beta")) %>% 
  filter(
    gamma %in% c(-0.2, 0, 0.2), marginal == 1,  beta %in% c(0,0.5)
  )%>% bind_cols(tibble(char = "db"))
bstPlotTib <- bind_rows(bstPlotTib, bstPlotTibDb)
#create factors and plot----
##themes----
#for jrssb: textSize <- 25
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
##factors----
bstPlotTib$k <- factor(bstPlotTib$k, levels = c(0,2,3,1),
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
bstPlotTib$char <- factor(bstPlotTib$char, 
                          levels = c("sb", "db"), 
                          labels = c(parse(text = TeX("$\\sigma^2_{sb}$")),
                                     parse(text = TeX("$\\sigma^2_{db}$"))
                          ))
#bstPlotTibDb <- bstPlotTib %>% filter(k== "db")
#bstPlotTib <- bstPlotTib %>% filter(k != "db")

rescale <- 100

#to have all values in the legend
dummy_data <- data.frame(
  m = c(NA, 0.2, NA, 0.4),  # Use NA to ensure these points are not visible in the plot
  varBst = c(0.1, NA, 0.3, NA),
  k = factor(c(0,1, 2, 3), levels = c(0,2,3,1),
             labels = c(parse(text = TeX("sb")), 
                        parse(text = TeX("cb(2)")), 
                        parse(text = TeX("cb(3)")), 
                        parse(text = TeX("db")))
  )
)

#Confidence intervals----
load(here("scripts/mean/bst/data/ciMeanTibRfix_Cheap"))
ciMeanTibRfixCheap <- ciMeanTibRfix
load(here("scripts/mean/bst/data/ciMeanTibRfix"))
#create tibbles----

#create factors----
ciMeanTibRfixCheap <- ciMeanTibRfixCheap %>% filter(
  beta %in% c(0, 0.5)
)
ciMeanTibRfix <- ciMeanTibRfix %>% filter(
  beta %in% c(0, 0.5)
)
ciMeanTibRfixCheap <- ciMeanTibRfixCheap %>% filter(
  gamma %in% c(-0.2,0, 0.2)
)
ciMeanTibRfix <- ciMeanTibRfix %>% filter(
  gamma %in% c(-0.2,0, 0.2)
)

ciMeanTibRfixCheap$k <- factor(ciMeanTibRfixCheap$k, levels = c(2,3,0,1),
                          labels = c(parse(text = TeX("cb(2)")), 
                                     parse(text = TeX("cb(3)")), 
                                     parse(text = TeX("sb")),
                                     parse(text = TeX("db"))
                          )
)
ciMeanTibRfix$k <- factor(ciMeanTibRfix$k, levels = c(2,3,0,1),
                               labels = c(parse(text = TeX("cb(2)")), 
                                          parse(text = TeX("cb(3)")), 
                                          parse(text = TeX("sb")),
                                          parse(text = TeX("db"))
                               )
)


ciMeanTibRfixCheap$marginal <- factor(ciMeanTibRfixCheap$marginal,
                                 levels = c(1,5),
                                 labels = c("GPD", "GEV")
)
ciMeanTibRfix$marginal <- factor(ciMeanTibRfix$marginal,
                                      levels = c(1,5),
                                      labels = c("GPD", "GEV")
)


ciMeanTibRfixCheap$beta <- factor(ciMeanTibRfixCheap$beta, 
                             levels = c(0, 0.25, 0.5, 0.75),
                             labels = c(parse(text = TeX("i.i.d.")), 
                                        parse(text = TeX("ARMAX $0.25$")), 
                                        parse(text = TeX("ARMAX $0.5$")), 
                                        parse(text = TeX("ARMAX $0.75$"))
                             )
)
ciMeanTibRfix$beta <- factor(ciMeanTibRfix$beta, 
                                  levels = c(0, 0.25, 0.5, 0.75),
                                  labels = c(parse(text = TeX("i.i.d.")), 
                                             parse(text = TeX("ARMAX $0.25$")), 
                                             parse(text = TeX("ARMAX $0.5$")), 
                                             parse(text = TeX("ARMAX $0.75$"))
                                  )
)


ciMeanTibRfixCheap$gamma <- factor(ciMeanTibRfixCheap$gamma, 
                              levels = c(0-0.2, -0.1, 0, 0.1, 0.2), 
                              labels = c(parse(text = TeX("$\\gamma = -0.2$")),
                                         parse(text = TeX("$\\gamma = -0.1$")),
                                         parse(text = TeX("$\\gamma = 0$")),
                                         parse(text = TeX("$\\gamma = 0.1$")),
                                         parse(text = TeX("$\\gamma = 0.2$"))
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
#create pdf with empirical coverage and relative average width 
# for different values of B for the CHEAP BOOTSTRAP
BVec <- c(1, 2, 5, 10, 20)
#start device for looping trough BVec
pdf(here("results/Cheap_Bootstrap_all_B_without_comparison.pdf"), 
    width = 11, height = 15)
for(loopB in BVec){
  ciCovPlot <- 
    ciMeanTibRfixCheap %>%   
    filter(k != "sb", marginal !="GEV", B == loopB) %>% 
    ggplot()+
    geom_line(aes(x = m, y = empCov, col = k),
              linewidth = 1.1)+
    geom_hline(yintercept = 0.95, col = "black", linetype = "longdash",
               linewidth = 1.1)+
    facet_grid(gamma~beta,  labeller = label_parsed)+
    labsPlot+
    labs( y = paste0("Empirical coverage"),
          title = paste0("Cheap Bootstrap CIs: B = ", loopB) )+
    #scale_y_continuous(
    #  breaks = c(0.91, 0.93, 0.95), limits = c(0.88, 0.98)
    #)+
    themePlot+
    scale_color_manual(values = ownPalette)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
    )
  #ciCovPlot
  
  
  ##now relative width Plot
  ciWidthDb <- ciMeanTibRfixCheap %>% filter(k == "db", B == loopB) %>% 
    select(-c(k,empCov)) %>% 
    rename(dbWidth = "avgWidth")
  ciWidthRel <- left_join(ciMeanTibRfixCheap %>% filter(B == loopB), ciWidthDb) %>% 
    filter(k %in% c("cb(2)", "cb(3)")) %>% mutate(relWidth = dbWidth/avgWidth)
  ciWidthRelPlot <- ciWidthRel %>% filter(marginal != "GEV") %>% 
    ggplot()+
    geom_line(aes(x = m, y = relWidth, col = k), linewidth = 1.1)+
    geom_hline(yintercept = 1, col = "#00BFC4", linetype = "longdash",
               linewidth = 1.1)+
    facet_grid(gamma~marginal+beta, scales = "free_y", labeller = label_parsed)+
    labsPlot+
    labs( y = paste0("Relative average width"),
          title = "Average relative width of confidence intervals",
    )+
    #scale_y_continuous(limits = c(1, 1.1), breaks = c(1,1.05, 1.1))+
    themePlot+
    scale_color_manual(values = ownPalette)+
    theme(strip.text.x = element_blank(), plot.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  #ciWidthRelPlot
  combCiPlot <- 
    ggarrange(ciCovPlot, ciWidthRelPlot, nrow = 2, align = "v", common.legend = T, 
              legend = "right",
              heights = c(1, 1.15))
  print(combCiPlot)
}
dev.off()





# as for B != 1 the performance of the Cheap BST is bad we only compare with 
# simple confidence intervals from our paper

# import the ci tibble for the simple CIs, create new col declaring if 
# the bstr_type of CI is "simple" or "cheap", filter only B == 1 in the Cheap tibble
# (as for B != 1 performance is bad), remove the col B, and bind_rows.
# Finally: plot with linetype = bstr_type

ciMeanTibRfix <- 
  ciMeanTibRfix %>% bind_cols(tibble(bstr_type = "simple"))

#start device for looping trough BVec
pdf(here("results/comparison_simple_cheap.pdf"), 
    width = 11, height = 15)
for (loopB in BVec){
  ciMeanTibRfixCheapB1 <- ciMeanTibRfixCheap %>% filter(B == loopB) %>% 
    select(-B) %>% bind_cols(bstr_type = "cheap")
  
  ciMeanTibRfix_both <- ciMeanTibRfix %>% bind_rows(ciMeanTibRfixCheapB1)
  
  ciCovPlot <- 
    ciMeanTibRfix_both %>%   
    filter(k != "sb", marginal !="GEV") %>% 
    ggplot()+
    geom_line(aes(x = m, y = empCov, col = k, linetype = bstr_type),
              linewidth = 1.1)+
    geom_hline(yintercept = 0.95, col = "black", linetype = "longdash",
               linewidth = 1.1)+
    facet_grid(gamma~beta,  labeller = label_parsed)+
    labsPlot+
    labs( y = paste0("Empirical coverage"),
          title = paste0("Comparison: Simple and Cheap (B= ", loopB, ") Bootstrap CIs")
    )+
    #scale_y_continuous(
    #  breaks = c(0.91, 0.93, 0.95), limits = c(0.88, 0.98)
    #)+
    themePlot+
    scale_color_manual(values = ownPalette)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
    )
  #ciCovPlot
  
  #rescale <- 10
  ciWidthPlot <-
    ciMeanTibRfix_both %>% 
    filter(k != "sb", marginal != "GEV") %>% 
    ggplot()+
    geom_line(aes(x = m, y = avgWidth, col = k, linetype = bstr_type),
              linewidth = 1.1)+
    facet_grid(gamma~beta+marginal, scales = "free_y", labeller = label_parsed)+
    labsPlot+
    labs( y = paste0("Average CI width"),
          title = "Average absolute width of confidence intervals",
    )+
    #scale_y_continuous(limits = c(1, 1.1), breaks = c(1,1.05, 1.1))+
    themePlot+
    scale_color_manual(values = ownPalette)+
    theme(strip.text.x = element_blank(), plot.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  #ciWidthPlot
  
  combCiPlot <- 
    ggarrange(ciCovPlot, ciWidthPlot, nrow = 2, align = "v", common.legend = T, 
              legend = "right",
              heights = c(1, 1.15))
  print(combCiPlot)
}
dev.off()
