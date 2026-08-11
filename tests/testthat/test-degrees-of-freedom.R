test_that("stage degrees of freedom use backend-native metadata", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  dat <- make_wine_test_data()

  fit <- tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  orm_df <- attr(stats::logLik(fit$fit_stage1), "df", exact = TRUE)
  vglm_rank <- methods::slot(fit$fit_stage2, "rank")

  expect_equal(.tsco_df(fit$fit_stage1), orm_df)
  expect_equal(.tsco_df(fit$fit_stage2), vglm_rank)
  expect_equal(
    attr(stats::logLik(fit), "df", exact = TRUE),
    orm_df + vglm_rank
  )

  lower_rank_vglm <- fit$fit_stage2
  methods::slot(lower_rank_vglm, "rank") <- as.integer(vglm_rank - 1L)

  expect_gt(length(stats::coef(lower_rank_vglm)), lower_rank_vglm@rank)
  expect_equal(.tsco_df(lower_rank_vglm), lower_rank_vglm@rank)
})


test_that("degree-of-freedom fallback excludes aliased coefficients", {
  dat <- data.frame(
    y = c(1, 3, 2, 5, 4, 7),
    x = seq_len(6)
  )
  dat$x_duplicate <- 2 * dat$x

  fit <- stats::lm(y ~ x + x_duplicate, data = dat)

  expect_lt(fit$rank, length(stats::coef(fit)))
  expect_equal(sum(is.finite(stats::coef(fit))), fit$rank)
  expect_equal(.tsco_df(fit), fit$rank)
})


test_that("failed ORM fits do not report model degrees of freedom", {
  fit <- structure(
    list(
      fail = TRUE,
      coefficients = c("Intercept" = 0, x = 1)
    ),
    class = c("orm", "rms")
  )

  expect_equal(.tsco_df(fit), NA_real_)
})
