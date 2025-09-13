library(here)
library(ggplot2)
library(tidyverse)
require(latex2exp)#for writing latex tau in the facets
textSize <- 20
source(here("backend/0plot_theme.R"))
plotPath <- file.path("Supplement", "plots")

dataPath <- here(file.path("backend", "data", "Supplement", "F.3."))
#-----------------------------------------------------------------------------
load(here(dataPath, "processedDataN5000.RData"))
varTibMleEst$marginal %>% unique()
varTibMleEst$tsMod %>% unique()
varTibMleEst$beta %>% unique()
#create tibbles----

varTibMleEstPl <- varTibMleEst %>% filter(beta %in% c(0,0.5)) 
#create factors and plots----
##factors----
varTibMleEstPl$k <- factor(varTibMleEstPl$k, levels = c(1,2,3,0),
                           labels = c(parse(text = TeX("db")), 
                                      parse(text = TeX("cb")), 
                                      parse(text = TeX("cb(3)")), 
                                      parse(text = TeX("sb"))
                           )
)
ownPalette <- #based on dark2
  c("cb" = "#F8766D",  
    "sb-cb"= "#F8766D",  
    "cb(2)" = "#F8766D",  
    "cb(3)" = "#7CAE00",  
    "db" = "#00BFC4",  
    "sb" = "#C77CFF")

varTibMleEstPl$marginal <- factor(varTibMleEstPl$marginal,
                                  levels = c(3,4),
                                  labels = c("Fréchet", "Pareto")
)

varTibMleEstPl$beta <- factor(varTibMleEstPl$beta, 
                              levels = c(0, 0.25, 0.5, 0.75),
                              labels = c(parse(text = TeX("i.i.d.")), 
                                         parse(text = TeX("ARMAX $0.25$")), 
                                         parse(text = TeX("ARMAX $0.5$")), 
                                         parse(text = TeX("ARMAX $0.75$"))
                              )
)

varTibMleEstPl$alpha <- factor(varTibMleEstPl$alpha, 
                               levels = c(0.5, 1, 1.5), 
                               labels = c(parse(text = TeX("$\\alpha_0 = 0.5$")),
                                          parse(text = TeX("$\\alpha_0 = 1$")),
                                          parse(text = TeX("$\\alpha_0 = 1.5$")))
)



##plots----
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
  col = "Estimator:")
rescale <- 100

#for jrssb
rescale <- 10
varTibMleEstPl$k <- factor(varTibMleEstPl$k,
                           levels = c("cb", "db", "cb(3)", "sb"),
                           labels = c("cb", "db", "cb(3)", "sb"))
msePlotMainJrssb <- 
  varTibMleEstPl %>% filter(marginal == "Pareto", k %in% c("db", "cb", "sb")) %>% 
  ggplot(aes(x = m, y = mse*rescale, col = k))+
  geom_line(linewidth = 1.1)+
  facet_grid(alpha~beta, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("MSE * ", rescale),
        title = "MSE of shape estimators")+
  themePlot+
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "right",
    plot.title = element_blank()
  )+
  scale_x_continuous(breaks = c(25,50,75),
                     limits = c(25,83)
  )+
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 2)
  )+
  scale_color_manual(values = ownPalette)

msePlotMainJrssb

if(TRUE){
  ggsave(msePlotMainJrssb, path = plotPath, 
         filename = "Figure F.8.pdf", device = "pdf",
         width = 13, height = 7, 
  )
}
