test_that("predict.tsco works with newdata and matches VGAM", {
  skip_if_not_installed("VGAM")

  data(pneumo, package = "VGAM")
  pneumo <- transform(pneumo, let = log(exposure.time))

  fit <- tsco(
    cbind(normal, mild, severe) ~ let,
    data = pneumo,
    levels = c("mild", "normal", "severe"),
    cutoff_level = "severe",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  nd <- data.frame(let = seq(min(pneumo$let), max(pneumo$let), length.out = 5))

  p_tsco <- predict(fit, newdata = nd, type = "prob")

  # Check structure
  expect_equal(nrow(p_tsco), nrow(nd))
  expect_equal(colnames(p_tsco), fit$levels)
  expect_equal(rowSums(p_tsco), rep(1, nrow(nd)), tolerance = 1e-8);

  # Compare to hand-calcuated predictvglm() calls
  p1 <- VGAM::predictvglm(fit$fit_stage1, newdata = nd, type = "response")
  p2 <- VGAM::predictvglm(fit$fit_stage2, newdata = nd, type = "response")

  p1 <- .tsco_align_prob(p1, fit$lower_levels)
  p2 <- .tsco_align_prob(p2, fit$collapsed_levels)

  p_manual <- matrix(
    0,
    nrow = nrow(nd),
    ncol = length(fit$levels),
    dimnames = list(NULL, fit$levels)
  )

  p_manual[, fit$lower_levels] <-
    p1[, fit$lower_levels, drop = FALSE] *
    as.numeric(p2[, fit$lower_collapsed_label])

  p_manual[, fit$upper_levels] <-
    p2[, fit$upper_levels, drop = FALSE]

  p_manual <- p_manual / rowSums(p_manual)

  expect_equal(
    p_tsco,
    p_manual,
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})
