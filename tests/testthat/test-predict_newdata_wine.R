test_that("predict.tsco works with newdata and matches VGAM", {
  dat <- make_wine_test_data()

  # This is the manuscript's headline PO|C|MR variant: conditional PO
  # below the cutoff, marginal (saturated) multinomial at/above it.
  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  nd <- dat$grouped[,c("temp", "contact")]

  p_tsco <- predict(fit, newdata = nd, type = "prob")

  # Check structure
  expect_equal(nrow(p_tsco), nrow(nd))
  expect_equal(colnames(p_tsco), fit$levels)
  expect_equal(rowSums(p_tsco), rep(1, nrow(nd)), tolerance = 1e-8, ignore_attr = TRUE)

  # Compare to hand-calculated predictvglm() calls
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

  expect_error(
    predict(fit, newdata = nd, backend_option = TRUE),
    regexp = "Additional prediction arguments are not supported"
  )
})
