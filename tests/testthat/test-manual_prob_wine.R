test_that("predict.tsco matches manual VGAM probability combination", {
  dat <- make_wine_test_data()

  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  nd <- expand.grid(
    temp = levels(dat$individual$temp),
    contact = levels(dat$individual$contact)
  )

  p_tsco <- predict(fit, newdata = nd, type = "prob")

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
