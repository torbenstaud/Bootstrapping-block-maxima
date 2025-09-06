facet_labels_gamma <- as_labeller(
  c(
    "-0.2" = "gamma == -0.2", 
    "-0.1" = "gamma == -0.1", 
    "0"    = "gamma == 0", 
    "0.1"  = "gamma == 0.1", 
    "0.2"  = "gamma == 0.2"
  ),
  label_parsed
)
facet_labels_beta <- 
  as_labeller(
    c("0" = "i.i.d.",
      "0.5" = "ARMAX 0.5"))


facet_label_targetvar <- as_labeller(
  c(
    "sb" = "sigma[sb]^{2}", 
    "db" = "sigma[db]^{2}"
  ),
  label_parsed
)