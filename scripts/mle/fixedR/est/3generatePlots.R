library(here)
library(Rcpp)
library(dplyr)
library(tidyr)
library(ggplot2)
require(latex2exp)#for writing latex tau in the facets
library(ggpubr)
#create QQ-Plots----
r <- 90
#vergleiche schätzer----
load(here(paste0("scripts/mle/fixedR/est/data/processedDataR",r,".RData")))

varTibMleEst <- varTibMleEst %>% 
  filter(marginal == 4, beta %in% c(0, 0.5)) #discard the Frechet marginals case
#create tibbles----
## absolute ----
varTibMleEstAbs <- varTibMleEst %>% 
  mutate(biasEst = biasEst^2) %>% rename(biasSquared = "biasEst") %>% 
  pivot_longer(
  cols = c("varEst", "biasSquared", "mse"), 
  names_to = "char", values_to = "val"
)

##relative mse----
dbTib <- varTibMleEst %>% filter(k == 1) %>% 
  select(-c(k,varEst, biasEst)) %>% 
  rename(dbMse = "mse")
  
varTibMleEstRel <- left_join(varTibMleEst, dbTib) %>% 
  filter(k != 1) %>% 
  mutate(ratioMse = dbMse / mse, relBias = biasEst^2/mse*100) %>% 
  select(-c(varEst, biasEst, mse, dbMse))

#create factors and plots----
##factors----

varTibMleEstAbs$k <- factor(varTibMleEstAbs$k, levels = c(1,2,3,0),
                            labels = c(parse(text = TeX("db")), 
                                       parse(text = TeX("cb(2)")), 
                                             parse(text = TeX("cb(3)")), 
                                       parse(text = TeX("sb"))
                                       )
                            )
varTibMleEstRel$k <- factor(varTibMleEstRel$k, levels = c(1,2,3,0),
                            labels = c(parse(text = TeX("db")), 
                                       parse(text = TeX("cb(2)")), 
                                       parse(text = TeX("cb(3)")), 
                                       parse(text = TeX("sb"))
                            )
                            )

varTibMleEstAbs$beta <- factor(varTibMleEstAbs$beta, 
                               levels = c(0, 0.25, 0.5, 0.75),
                               labels = c(parse(text = TeX("i.i.d.")), 
                                          parse(text = TeX("ARMAX $0.25$")), 
                                          parse(text = TeX("ARMAX $0.5$")), 
                                          parse(text = TeX("ARMAX $0.75$"))
                               )
                               )
varTibMleEstRel$beta <- factor(varTibMleEstRel$beta, 
                               levels = c(0, 0.25, 0.5, 0.75),
                               labels = c(parse(text = TeX("i.i.d.")), 
                                          parse(text = TeX("ARMAX $0.25$")), 
                                          parse(text = TeX("ARMAX $0.5$")), 
                                          parse(text = TeX("ARMAX $0.75$"))
                               )
                               )
varTibMleEstAbs$alpha <- factor(varTibMleEstAbs$alpha, 
                                levels = c(0.5, 1, 1.5), 
                                labels = c(parse(text = TeX("$\\alpha = 0.5$")),
                                           parse(text = TeX("$\\alpha = 1$")),
                                           parse(text = TeX("$\\alpha = 1.5$")))
                                )

varTibMleEstRel$alpha <- factor(varTibMleEstRel$alpha, 
                                levels = c(0.5, 1, 1.5), 
                                labels = c(parse(text = TeX("$\\alpha_0 = 0.5$")),
                                           parse(text = TeX("$\\alpha_0 = 1$")),
                                           parse(text = TeX("$\\alpha_0 = 1.5$")))
)


##plots----
textSize <- 12
themePlot <- theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2),
                   strip.background = element_rect(color = "black", 
                                                   fill = "lightgrey", size = 0.2),
                   axis.title.x = element_text(size = textSize),
                   axis.title.y = element_text(size = textSize),
                   axis.text.y =element_text(size=textSize), 
                   axis.text.x =element_text(size=textSize, angle = 90, 
                                             hjust = 1, vjust = 0.5),
                   strip.text.x = element_text(size = textSize),
                   strip.text.y = element_text(size = textSize),
                   plot.title = element_text(hjust = 0.5, size = textSize, 
                                             face = "bold"), 
                   #panel.background = element_rect(rgb(0.95, 0.95, 0.95, alpha = 1)),
                   legend.position = "right",
                   legend.title = element_text(size = textSize),
                   legend.text = element_text(size = textSize))

labsPlot <- labs(title = "Fréchet Shape estimation (r = 90 fixed)",
                 x = "Effective sample size m",
                 col = "Estimator")

plotRelMse <- 
  varTibMleEstRel %>% ggplot(aes(x = m, y = ratioMse, col = k))+
  geom_line()+
  geom_hline(yintercept = 1, col = "black", linetype = "dashed")+
  facet_grid(alpha~beta, labeller = label_parsed)+ 
  labsPlot+
  labs( y = "Relative MSE")+
  themePlot+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
 
plotRelBias <- 
varTibMleEstRel %>% 
  ggplot(aes(x = m, y = relBias, col = k))+
  geom_line()+
  facet_grid(alpha~beta, labeller = label_parsed)+
  labs(y = "Bias^2/MSE * 100")+
  labsPlot +
  themePlot+
  theme(strip.background.x = element_blank(), 
        strip.text.x = element_blank(),
        plot.title =  element_blank())
  #scale_x_continuous(breaks = c(50,60, 7080 , 100))+
  #ylim(0,0.8)+
  #scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6))+

plotRelMseBias <- ggarrange(plotRelMse, plotRelBias, align = "v", ncol = 1,
                            common.legend = T, legend = "bottom")
plotRelMseBias
if(FALSE){
  ggsave(plotRelMseBias, filename = "plotFreFixRRelMseBias.pdf",
         device = "pdf", path = here("results/"),
         width = 6, height = 10)
}

