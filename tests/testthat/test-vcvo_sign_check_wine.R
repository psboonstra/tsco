test_that("vcov.tsco applies documented sign correction for all stage combinations", {
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

    V2_raw <- as.matrix(stats::vcov(fit$fit_stage2))

    s2 <- .tsco_sign_vector(
      fit$fit_stage2,
      kind = fit$stage2,
      po.reverse = fit$po.reverse
    )

    if (!is.null(rownames(V2_raw))) {
      s2 <- s2[rownames(V2_raw)]
    }

    V2_expected <- V2_raw * outer(s2, s2)

    expect_equal(
      vcov(fit, stage = "stage2"),
      V2_expected,
      tolerance = 1e-10
    )
  }
})
