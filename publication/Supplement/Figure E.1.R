# --------------------------------------------------------------------
# Script Name: meanAsyVar.R
# Purpose: This scripts produces the relative asymptotic variance plot
#   for the mean estimator
# Author: Torben Staud
# Date: 2024-08
# --------------------------------------------------------------------

# Clear the environment and console
rm(list=ls())  # Clear all objects from the current R environment
cat("\014")    # Clear the console

# ==== SECTION 1: Load required libraries ======================================
# Description: Load all the R packages required for the script.
#_______________________________________________________________________________
library(Rcpp)
library(here)
library(tidyverse)
library(pracma)# integration
library(tidyverse) #tidyverse package suite for data manipulation

section <- "Supplement"


# ==== SECTION 2: main script ===========================================
# Description: code for the formulas and plot
#______________________________________________________________________________
textSize <- 20
source(here("backend/0plot_theme.R"))
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

# ==== SECTION 3: main script ===========================================
# Description: code for the formulas and plot
#______________________________________________________________________________

##disjoints blocks----
var.db <- function(gam){
  ifelse(gam == 0,
         pi^2/6,
         (gamma(1-2*gam)-gamma(1-gam)^2)/gam^2)
}

## sliding blocks----
#auxiliary function representing alpha_(2*gamma)(x)
h.f <- function(gam,x) {
  b <- (2)*gam
  if (b != -1) {
    v <- (1-(1-x)^(b +1))/(x*(b+1))
  }
  else {
    v <- 
      -log(1-x)/x
    
  }
  return(v)
}


C.int <- function(gam){
  if(sum(gam >= 1/2) > 0) {
    stop("gamma has to be smaller than 1/4")
  }
  integrand <- function(x) {
    v <- 
      (h.f(gam,x)-1) * ((1-x)*x)^(-gam-1)
    
    return(v)
  }
  return(
    quadgr(integrand, 0, 1/2)$value
  )
}


var.sl <- function(gam_vec) {
  val <- numeric(length(gam_vec))
  for ( i in seq(1:length(val))) {
    gam <- gam_vec[i]
    if (gam_vec[i] > 0) {
      val[i] <- 
        -2*gamma(1-2*gam)/gam*C.int(gam)
    }
    else if(gam_vec[i] < 0) {
      val[i] <-
        4*gamma(-2*gam)*C.int(gam)
    }
    else {
      val[i] <-
        4*(log(4)-1)
    }
  }
  return(
    val
  )
} 

##plots
xx <- seq(-0.75, 1/2-10^(-3), length.out = 100)
var_db_dat <- var.db(xx)
var_sl_dat <- var.sl(xx)
quot <- var_sl_dat/var_db_dat
plt_data <- tibble(type = rep(c("db","sl"), each = length(xx)), 
                   gamma = c(xx,xx),
                   sigma = c(var_db_dat,var_sl_dat)
)
plt_quot_data <- tibble(gamma = xx, ratio = quot)

p_var <- plt_data %>% ggplot(aes(x = gamma, y = sigma, col = as.factor(type)))+
  geom_line()
p_var
p_quot <- plt_quot_data %>% ggplot(aes(x = gamma, y = ratio))+
  labs(title = "Mean estimation",
       x = expression(gamma),
       y = "Ratio")+
  geom_line()+
  geom_hline(
    yintercept = 1, col = "red", linetype = "dashed"
  )+
  themePlot+
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  scale_x_continuous(limits = c(-0.5, 0.5))
#scale_y_continuous(limits = c(1, 1.17))
p_quot
if(T){
  ggsave(p_quot, filename = "Figure E.1.pdf", device = "pdf", 
         path  = file.path(here(section), "plots"),
         height = 4, width = 8)
}


