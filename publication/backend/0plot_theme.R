#note: need to define the variable "textSize" before sourcing this theme file
##often values: 20 - 25 are good

themePlot <- theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2),
                   strip.background = element_rect(color = "black", 
                                                   fill = "lightgrey", size = 0.2),
                   axis.title.x = element_text(size = textSize),
                   axis.title.y = element_text(size = textSize),
                   axis.text.y =element_text(size=textSize), 
                   axis.text.x =element_text(size=textSize, angle = 90),
                   strip.text.x = element_text(size = textSize),
                   strip.text.y = element_text(size = textSize),
                   plot.title = element_text(hjust = 0.5, size = textSize, 
                                             face = "bold"), 
                   #panel.background = element_rect(rgb(0.95, 0.95, 0.95, alpha = 1)),
                   legend.position = "right",
                   legend.title = element_text(size = textSize),
                   legend.text = element_text(size = textSize))
ownPalette <- #based on dark2
  c("cb(2)" = "#F8766D",
    "cb" = "#F8766D",
    "cb(3)" = "#7CAE00",  
    "db" = "#00BFC4",  
    "sb" = "#C77CFF",
    "sb-cb" = "#F8766D")

facet_labels_param <- as_labeller(c(
  "location" = "mu",
  "shape" = "gamma",
  "scale" = "sigma",
  "return_level(100,365)" = "RL"
), label_parsed)
color_db <- "#00BFC4"
color_cb <- "#F8766D"
