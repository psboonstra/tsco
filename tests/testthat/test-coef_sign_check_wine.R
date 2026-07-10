test_that("coef.tsco applies documented sign correction for all stage combinations", {
  dat <- make_wine_test_data()

  combos <- expand.grid(
    stage1 = c("po", "multinomial"),
    stage2 = c("po", "multinomial"),
    stringsAsFactors = FALSE
  )

  for (ii in seq_len(nrow(combos))) {
    fit <- tsco(
      dat$grouped_formula,
      data = dat$grouped,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = combos$stage1[ii],
      stage2 = combos$stage2[ii],
      warn_degenerate = FALSE
    )

    b1_raw <- stats::coef(fit$fit_stage1)
    b2_raw <- stats::coef(fit$fit_stage2)

    s1 <- .tsco_sign_vector(
      fit$fit_stage1,
      kind = fit$stage1,
      po.reverse = fit$po.reverse
    )

    s2 <- .tsco_sign_vector(
      fit$fit_stage2,
      kind = fit$stage2,
      po.reverse = fit$po.reverse
    )

    expect_equal(coef(fit, stage = "stage1"), b1_raw * s1)
    expect_equal(coef(fit, stage = "stage2"), b2_raw * s2)
  }
})
