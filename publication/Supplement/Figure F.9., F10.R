library(here)
library(ggplot2)
library(tidyverse)
require(latex2exp)#for writing latex tau in the facets
library(ggpubr)

textSize <- 20
source(here("backend/0plot_theme.R"))
plotPath <- file.path("Supplement", "plots")

dataPath <- here(file.path("backend", "data", "Supplement", "F.3."))
#-----------------------------------------------------------------------------
load(file.path(dataPath, "trueVarsTib"))
load(file.path(dataPath, "bstFreVarTib"))
#create tibbles----

bstFreVarTibPl <- bstFreVarTib %>% 
  left_join(trueVarsTib) %>% 
  filter( beta %in% c(0, 0.5), k != 1) %>% 
  bind_cols(tibble(char = "sb"))


## for disjoint blocks----
load(file.path(dataPath, "trueVarsTibDb"))#db true vars
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
  c("sb-cb" = "#F8766D",
    "cb" = "#F8766D",
    "cb(2)" = "#F8766D",  
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
labsPlot <- labs(
  x = "Effective sample size m",
  col = "Bootstrap:")
rescale <- 100
#now for the main paper: no Fréchet marginals-----
bstFreVarTibPl$k <- factor(bstFreVarTibPl$k, 
                           levels = c("db", "sb", "cb(2)", "cb(3)"),
                           labels = c("db", "sb", "cb", "cb(3)")
                           
)
bstVarPlotMain <- 
  bstFreVarTibPl %>% filter(marginal == "Pareto", 
                            k %in% c("db", "cb")) %>% 
  mutate(asVar = m*varBst, asTrueVar = m*trueVar) %>% 
  filter(alpha != 2) %>% 
  ggplot()+
  geom_line(aes(x = m, y = varBst*rescale, col = k),
            linewidth = 1.1)+
  geom_line(aes(x = m, y = trueVar*rescale), col = "black", linetype = "dashed",
            linewidth = 1.1)+
  facet_grid(alpha~char+beta, scales = "free_y", labeller = label_parsed)+
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
#for jrssb
if(TRUE){
  ggsave(bstVarPlotMain, path = plotPath, 
         filename = "Figure F.9.pdf", device = "pdf",
         width = 13, height = 7, 
  )
}


#Confidence Intervals----
load(file.path(dataPath, "CiFreNfix"))
ciFreTibNfix
##sizes----
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
##create factors----
ciFreTibNFixPlot <- ciFreTibNfix %>% filter(
  beta %in% c(0, 0.5), 
  m >=25, m<=83
)




ciFreTibNFixPlot$k <- factor(ciFreTibNFixPlot$k, levels = c(0,2,3,1),
                             labels = c(parse(text = TeX("sb")),
                                        parse(text = TeX("cb")), 
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

#now for the main paper, where we omit the Fréchet marginal case
ciFreTibNFixPlotMain <- ciFreTibNFixPlot %>% filter(marginal == "Pareto")
ciCovPlotMain <- 
  ciFreTibNFixPlotMain %>%  filter( marginal == "Pareto") %>% 
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
  filter(k %in% c("cb", "cb(3)")) %>% mutate(relWidth = avgWidth/dbWidth)

#for jrssb
ciFreTibNFixPlotMain$k <- 
  factor(ciFreTibNFixPlotMain$k, 
         levels = c("sb", "db", "cb", "cb(3)"),
         labels = c("sb", "db", "sb-cb", "cb(3)"))
ciCovPlotMainJrssb <- 
  ciFreTibNFixPlotMain %>%  filter( marginal == "Pareto") %>% 
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = empCov, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 0.95, col = "black", linetype = "dashed",
             linewidth = 1.1)+
  facet_grid(alpha~beta, scales = "free_y", labeller = label_parsed)+
  themePlot+
  labsPlot+
  labs(col = "Confidence Interval:")+
  labs( y = paste0("Empirical coverage"),
        title = "Confidence intervals for the shape")+
  scale_color_manual(values = ownPalette)+
  theme(strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_blank()
  )+
  scale_y_continuous(
    limits = c(0.85, 1)
  )+
  scale_x_continuous(breaks = c(25, 50, 75))
ciCovPlotMainJrssb


ciWidthRelPlotMainJrssb <- ciWidthRelMain %>% filter( ) %>% 
  ggplot()+
  geom_line(aes(x = m, y = relWidth, col = k),
            linewidth = 1.1)+
  geom_hline(yintercept = 1, col = "#00BFC4",
             linewidth = 1.1)+
  facet_grid(alpha~beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs(col = "Confidence Interval:")+
  labs( y = paste0("Relative average width"),
        title = "Average relative width of confidence intervals",
  )+
  #scale_y_continuous(limits = c(0.95, 1.25), breaks = c(1,1.1,1.2))+
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

if(TRUE){
  ggsave(combCiPlotMainJrssb, path = plotPath, 
         filename = "Figure F.10.pdf", device = "pdf",
         width = 13, height = 7
  )
}
