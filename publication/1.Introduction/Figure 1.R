library(here)
library(tidyverse)
textSize <- 20
source(here("backend/0plot_theme.R"))
section <- "1.Introduction"
#analyze bstr/est errors by means of histograms
gamma <- -0.2
beta <- 0.5
m <- 80
r <- 365
N <- 1000
B <- 1000
nameLoad <- 
  sprintf("array_bstEst_gam%.1f_beta%.1f_m%g_r%g", gamma, beta, m, r)
load(file = file.path("backend/data/1", nameLoad)) #150mb file - 
load("backend/data/1/truth_array.dat")
truth <- array_truth[1, 2, 4]

array_estErrs <- 
  array_bsts[,1,4,] - mean(array_bsts[,1,4,])


#bootstrap errors
array_bstErrs <- 
  apply(array_bsts[,-c(B+2),4,], c(1,3), function(bsts) bsts[-1] - bsts[1]) %>% 
  aperm(c(2,3,1)) %>% 
  apply(c(2), function(bsts) c(bsts))


tib_full_errs <- 
  bind_rows(
    tibble(x = array_estErrs[,1], method = "sb", char = "est_err"),
    tibble(x = array_estErrs[,2], method = "db", char = "est_err"),
    tibble(x = array_estErrs[,3], method = "cb", char = "est_err"),
    tibble(x = array_bstErrs[seq(1,N),1], method = "sb", char = "bst_err") %>% 
      filter(x <= max(array_estErrs[,1]), x >= min(array_estErrs[,1])),
    tibble(x = array_bstErrs[seq(1,N),2], method = "db", char = "bst_err") %>% 
      filter(x <= max(array_estErrs[,2]), x >= min(array_estErrs[,2])),
    tibble(x = array_bstErrs[seq(1,N),3], method = "cb", char = "bst_err") %>% 
      filter(x <= max(array_estErrs[,3]), x >= min(array_estErrs[,3]))
  )

tib_full_errs <- tib_full_errs %>% rename(value = x)
tib_full_errs$method <- factor(tib_full_errs$method,
                               levels = c("sb", "cb", "db"))

facet_labels_method <- as_labeller(
  c(
    "db" = "plain('err')^'[d]' ~ 'vs' ~ plain('err')^'[d],*'",
    "sb" = "plain('err')^'[s]' ~ 'vs' ~ plain('err')^'[s],*'",
    "cb" = "plain('err')^'[s]' ~ 'vs' ~ plain('err')^'[c],*'"
  ),
  label_parsed
)

## for the intro:
facet_labels_method2 <- as_labeller(
  c(
    "db" = "Disjoint",
    "sb" = "Naive Sliding",
    "cb" = "Circular"
  )
)
plotRlTs_histsIntro <-
  tib_full_errs %>% 
  filter(method != "db") %>% 
  ggplot(aes(x = value, y= after_stat(density)))+
  geom_histogram(aes(x = value, fill = char), 
                 alpha = 0.4, position = "identity", bins = 20)+
  geom_density(aes(x = value, linetype = char), linewidth = 1)+
  facet_grid(cols = vars(method), labeller = labeller(method = facet_labels_method2))+
  labs(
    y = "Density",
    x = "Centered Errors",
    fill = "Error Type:", linetype = "Error Type:"
  )+
  scale_fill_manual(labels = c("Bootstrap", "Estimation"), values = c("#F8766D", "#00BFC4"))+
  scale_linetype_manual(labels = c("Bootstrap", "Estimation"), values = c("solid", "dashed"))+
  themePlot+
  theme(
    #axis.title.x = element_blank()
    legend.position = "right"
  )

plotRlTs_histsIntro
if(T){
  ggsave(plotRlTs_histsIntro, 
         filename = sprintf("Figure 1.pdf"),
         device = "pdf", path = file.path(here(section), "plots"),
         width = 13, height = 4)
}
