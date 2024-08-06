library(here)
library(ggplot2)
library(tidyverse)
require(latex2exp)#for writing latex tau in the facets
library(ggpubr)
#-----------------------------------------------------------------------------
load(here("scripts/mle/fixedR/bst/generateTrueVars/data/trueVarsTib"))
load(here("scripts/mle/fixedR/bst/data/bstFreVarTib"))
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
  #filter(k != "sb") %>% 
  ggplot()+
  geom_line(aes(x = m, y = varBst*rescale, col = k))+
  geom_line(aes(x = m, y = trueVar*rescale), col = "black", linetype = "dashed")+
  facet_grid(alpha~beta+marginal, scales = "free_y", labeller = label_parsed)+
  labsPlot+
  labs( y = paste0("Variance * ",rescale),
        title = "Bootstrapping the variance of the shape estimator (r = 90 fixed)")+
  themePlot
bstVarPlot
if(FALSE){
  ggsave(bstVarPlot, path = here("results/"), 
         filename = "plotFreFixRVarEst.pdf", device = "pdf",
         width = 10, height = 8, 
  )
}




