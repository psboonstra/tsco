
check_tsco <- function(probabilities_ctrl, probabilities_trt, digits = 2) {
  
  n_levels <- length(probabilities_ctrl)
  ordered_logits_ctrl <- qlogis(rev(cumsum(rev(probabilities_ctrl[-n_levels]) / (1 - probabilities_ctrl[n_levels])))[-1])
  ordered_logits_trt <- qlogis(rev(cumsum(rev(probabilities_trt[-n_levels]) / (1 - probabilities_trt[n_levels])))[-1])
  
  # This is the log-or for the last category
  foo2 = qlogis(probabilities_ctrl[n_levels]) - qlogis(probabilities_trt[n_levels])
  
  # check for one distinct value of ordered_logits_ctrl - ordered_logits_trt:
  foo = ordered_logits_ctrl - ordered_logits_trt
  as.character(ifelse(near(max(foo), min(foo)), paste0(round(mean(foo), digits),",", round(foo2, digits), collapse = ","), "No"))
  
}
