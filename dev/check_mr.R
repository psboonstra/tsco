# Identify odds ratios from control probabilities and treatment probabilities:
check_mr <- function(probabilities_ctrl, probabilities_trt) {
  stopifnot(all(probabilities_ctrl >= 0) && abs(sum(probabilities_ctrl) - 1) < .Machine$double.eps^0.5)
  stopifnot(all(probabilities_trt >= 0) && abs(sum(probabilities_trt) - 1) < .Machine$double.eps^0.5)
  stopifnot(length(probabilities_ctrl) == length(probabilities_trt))
  
  logits_ctrl <- log(probabilities_ctrl) - log(probabilities_ctrl)[1] 
  logits_trt <- log(probabilities_trt) - log(probabilities_trt)[1] 
  
  foo <- logits_ctrl - logits_trt
  foo
}