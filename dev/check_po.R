check_po <- function(probabilities_ctrl, probabilities_trt, digits = 2) {
  
  n_levels <- length(probabilities_ctrl)
  ordered_logits_ctrl <- qlogis(rev(cumsum(rev(probabilities_ctrl)))[-1])
  ordered_logits_trt <- qlogis(rev(cumsum(rev(probabilities_trt)))[-1])
  # check for one distinct value of ordered_logits_ctrl - ordered_logits_trt:
  foo = ordered_logits_ctrl - ordered_logits_trt
  as.character(ifelse(near(max(foo), min(foo)), round(mean(foo), digits), "No"))
}
