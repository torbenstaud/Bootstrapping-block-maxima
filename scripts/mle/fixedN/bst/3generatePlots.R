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
  filter( beta %in% c(0, 0.5), k != 1)



#create factors and plots----
##factors----
bstFreVarTibPl$k <- factor(bstFreVarTibPl$k, levels = c(1,2,3,0),
                           labels = c(parse(text = TeX("db")), 
                                      parse(text = TeX("cb(2)")), 
                                      parse(text = TeX("cb(3)")), 
                                      parse(text = TeX("sb"))
                           )
)
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

##plots----
# textSize = 15
labsPlot <- labs(
  x = "Effective sample size m",
  col = "Bootstrap")
rescale <- 100

bstVarPlot <- 
bstFreVarTibPl %>% mutate(asVar = m*varBst, asTrueVar = m*trueVar) %>% 
  filter(alpha != 2) %>% 
  ggplot()+
  geom_line(aes(x = m, y = varBst*rescale, col = k))+
  geom_line(aes(x = m, y = trueVar*rescale), col = "black", linetype = "dashed")+
  facet_grid(alpha~beta+marginal, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Variance * ",rescale),
        title = "Bootstrapping the variance of the shape estimator (n = 1,000 fixed)")+
  themePlot
bstVarPlot
if(FALSE){
  ggsave(bstVarPlot, path = here("results/"), 
         filename = "plotFreFixNVarEst.pdf", device = "pdf",
         width = 10, height = 8, 
  )
}





#Confidence Intervals----
load(here("scripts/mle/fixedN/bst/data/CiFreNfix"))
ciFreTibNfix
#create tibbles----

#create factors----
ciFreTibNFixPlot <- ciFreTibNfix %>% filter(
  beta %in% c(0, 0.5)
)




ciFreTibNFixPlot$k <- factor(ciFreTibNFixPlot$k, levels = c(2,3,0,1),
                             labels = c(parse(text = TeX("cb(2)")), 
                                        parse(text = TeX("cb(3)")), 
                                        parse(text = TeX("sb")),
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
  ciFreTibNFixPlot %>%   filter( m != 31) %>% 
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = empCov, col = k))+
  geom_hline(yintercept = 0.95, col = "black", linetype = "dashed")+
  facet_grid(alpha~beta+marginal, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Average empirical coverage"),
        title = "Confidence intervals for the shape")+
  themePlot+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
ciCovPlot

#rescale <- 10
ciWidthPlot <-
  ciFreTibNFixPlot %>%  
  filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = avgWidth, col = k))+
  facet_grid(alpha~beta+marginal, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Width"),
        title = "Average width of confidence intervals")+
  themePlot
ciWidthPlot

##now relative width Plot
ciWidthDb <- ciFreTibNFixPlot %>% filter(k == "db") %>% select(-c(k,empCov)) %>% 
  rename(dbWidth = "avgWidth")
ciWidthRel <- left_join(ciFreTibNFixPlot, ciWidthDb) %>% 
  filter(k %in% c("cb(2)", "cb(3)")) %>% mutate(relWidth = dbWidth/avgWidth)
ciWidthRelPlot <- ciWidthRel %>% filter( m != 31) %>% 
  ggplot()+
  geom_line(aes(x = m, y = relWidth, col = k))+
  geom_hline(yintercept = 1, col = "black", linetype = "dashed")+
  facet_grid(alpha~beta+marginal, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Relative average width"),
        title = "Average relative width of confidence intervals",
        )+
  scale_y_continuous(limits = c(0.95, 1.25), breaks = c(1,1.1,1.2))+
  themePlot+
  theme(strip.text.x = element_blank(), plot.title = element_blank())

ciWidthRelPlot

combCiPlot <- 
  ggarrange(ciCovPlot, ciWidthRelPlot, nrow = 2, align = "v", common.legend = T, 
            legend = "bottom")
combCiPlot

if(FALSE){
  ggsave(combCiPlot, path = here("results/"), 
         filename = "plotFreFixNCi.pdf", device = "pdf",
         width = 10, height = 10, 
  )
}
