
check_cppo <- function(probabilities_ctrl, probabilities_trt, digits = 2) {
  
  n_levels <- length(probabilities_ctrl)
  ordered_logits_ctrl <- qlogis(rev(cumsum(rev(probabilities_ctrl)))[-1])
  ordered_logits_trt <- qlogis(rev(cumsum(rev(probabilities_trt)))[-1])
  
  # check for one distinct value of ordered_logits_ctrl - ordered_logits_trt:
  foo = (ordered_logits_ctrl - ordered_logits_trt)[-(n_levels-1)]
  foo2 = (ordered_logits_ctrl - ordered_logits_trt)[(n_levels-1)] - mean(foo)
  as.character(ifelse(near(max(foo), min(foo)),  paste0(round(mean(foo), digits),",", round(foo2, digits), collapse = ","), "No"))
}
