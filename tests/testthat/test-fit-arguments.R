test_that("stage-specific fitter arguments are routed and reused", {
  dat <- make_wine_test_data()
  stage1_args <- list(x = TRUE, y = TRUE)
  stage2_args <- list(x.arg = FALSE, y.arg = FALSE)

  fit <- tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    stage1_args = stage1_args,
    stage2_args = stage2_args,
    warn_degenerate = FALSE
  )

  expect_identical(fit$stage1_args, stage1_args)
  expect_identical(fit$stage2_args, stage2_args)
  expect_equal(nrow(fit$fit_stage1$x), fit$n_stage1_obs)
  expect_equal(length(fit$fit_stage1$y), fit$n_stage1_obs)
  expect_equal(dim(fit$fit_stage2@x), c(0L, 0L))
  expect_equal(dim(fit$fit_stage2@y), c(0L, 0L))

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  expect_equal(sum(!is.finite(lrt$Chisq)), 0L)
})


test_that("stage-specific fitter arguments are validated", {
  dat <- make_wine_test_data()

  expect_snapshot(error = TRUE, {
    tsco(
      dat$individual_formula,
      data = dat$individual,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1_args = TRUE
    )
  })

  expect_snapshot(error = TRUE, {
    tsco(
      dat$individual_formula,
      data = dat$individual,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1_args = list(TRUE)
    )
  })

  expect_snapshot(error = TRUE, {
    tsco(
      dat$individual_formula,
      data = dat$individual,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage2_args = structure(
        list(TRUE, FALSE),
        names = c("x.arg", "x.arg")
      )
    )
  })

  expect_snapshot(error = TRUE, {
    tsco(
      dat$individual_formula,
      data = dat$individual,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage2_args = list(data = dat$individual, family = "logistic")
    )
  })
})
