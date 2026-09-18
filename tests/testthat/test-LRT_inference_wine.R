test_that("negative LRT statistics are cleaned or rejected", {
  tiny <- .tsco_clean_lrt(-1e-9, "x")
  material <- .tsco_clean_lrt(-1e-4, "x")

  expect_equal(tiny$value, 0)
  expect_identical(tiny$problem, NA_character_)
  expect_identical(material$value, NA_real_)
  expect_match(material$problem, "reduced model had a larger log-likelihood")
})


test_that("failed reduced LRT fits retain diagnostic information", {
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

  # Deliberately make reduced stage-1 refits fail. This modifies only the
  # stored refit data; the already fitted full model remains valid.
  fit$data_stage1$.tsco_y_stage1 <- NULL

  expect_warning(
    lrt <- summary(fit, joint_test = "LRT")$joint_tests,
    "Joint likelihood-ratio test problems"
  )

  # Both formula terms should have failed because every reduced stage-1 model
  # requires the removed response.
  expect_true(all(is.na(lrt$Chisq)))
  expect_true(all(is.na(lrt$df)))
  expect_true(all(is.na(lrt[["Pr(>Chisq)"]])))

  errors <- attr(lrt, "errors")

  expect_false(is.null(errors))
  expect_identical(names(errors), c("temp", "contact"))
  expect_length(errors, 2L)
  expect_true(all(nzchar(errors)))
})



test_that("LRT in summary matches manual reduced-model refit", {
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

  s <- summary(fit, joint_test = "LRT")

  expect_true(all(
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in%
      colnames(s$coef_stage1)
  ))

  expect_true(all(
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in%
      colnames(s$coef_stage2)
  ))

  expect_true(all(c("temp", "contact") %in% rownames(s$joint_tests)))
  expect_true(all(s$joint_tests[c("temp", "contact"), "df"] == 2))
  expect_true(all(is.finite(s$joint_tests[c("temp", "contact"), "Chisq"])))
  expect_true(all(is.finite(s$joint_tests[c("temp", "contact"), "Pr(>Chisq)"])))

  expect_true(all(is.finite(s$joint_tests$Chisq)))
  expect_true(all(s$joint_tests$Chisq >= -1e-8))
  expect_true(all(is.finite(s$joint_tests[["Pr(>Chisq)"]])))

  term <- attr(stats::terms(fit$stage1_formula), "term.labels")[2]

  f1_reduced <- .tsco_drop_term_formula(fit$stage1_formula, term)
  f2_reduced <- .tsco_drop_term_formula(fit$stage2_formula, term)

  fit1_reduced <- .tsco_fit_stage(
    formula = f1_reduced,
    kind = fit$stage1,
    engine = fit$stage1_engine,
    data = fit$data_stage1,
    weights = fit$weights_stage1,
    stage_args = fit$stage1_args
  )

  fit2_reduced <- .tsco_fit_stage(
    formula = f2_reduced,
    kind = fit$stage2,
    engine = fit$stage2_engine,
    data = fit$data_stage2,
    weights = fit$weights_stage2,
    stage_args = fit$stage2_args
  )

  ll_full <- .tsco_stage_loglik(fit$fit_stage1, fit$stage1_engine) +
    .tsco_stage_loglik(fit$fit_stage2, fit$stage2_engine)

  ll_reduced <- .tsco_stage_loglik(fit1_reduced, fit$stage1_engine) +
    .tsco_stage_loglik(fit2_reduced, fit$stage2_engine)

  lrt_manual <- 2 * (ll_full - ll_reduced)

  df_manual <- (
    .tsco_df(fit$fit_stage1) - .tsco_df(fit1_reduced)
  ) + (
    .tsco_df(fit$fit_stage2) - .tsco_df(fit2_reduced)
  )

  expect_equal(
    s$joint_tests[term, "Chisq"],
    lrt_manual,
    tolerance = 1e-8
  )

  expect_equal(
    s$joint_tests[term, "df"],
    df_manual
  )
})


test_that("mixed-backend LRT matches independent manual refits", {
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

  expect_identical(fit$stage1_engine, "orm")
  expect_identical(fit$stage2_engine, "vglm")

  term <- "contact"
  f1_reduced <- stats::update(fit$stage1_formula, paste(". ~ . -", term))
  f2_reduced <- stats::update(fit$stage2_formula, paste(". ~ . -", term))

  fit1_reduced <- rms::orm(
    f1_reduced,
    data = fit$data_stage1,
    family = "logistic"
  )
  fit2_reduced <- VGAM::vglm(
    f2_reduced,
    data = fit$data_stage2,
    family = VGAM::multinomial(refLevel = 1)
  )

  # Deliberately uses the backends' *native* log-likelihoods. For
  # individual-level data those coincide with the package's constant-free
  # scale (every row is one draw, so the multinomial constants are zero), and
  # this test is the independent check that they do. The grouped analogue
  # below covers the case where the two scales differ.
  ll_full <- as.numeric(stats::logLik(fit$fit_stage1)) +
    as.numeric(stats::logLik(fit$fit_stage2))
  ll_reduced <- as.numeric(stats::logLik(fit1_reduced)) +
    as.numeric(stats::logLik(fit2_reduced))

  lrt_manual <- 2 * (ll_full - ll_reduced)
  df_manual <-
    length(stats::coef(fit$fit_stage1)) - length(stats::coef(fit1_reduced)) +
    length(stats::coef(fit$fit_stage2)) - length(stats::coef(fit2_reduced))

  joint_tests <- summary(fit, joint_test = "LRT")$joint_tests

  expect_equal(
    joint_tests[term, "Chisq"],
    lrt_manual,
    tolerance = 1e-8
  )
  expect_equal(joint_tests[term, "df"], df_manual)
})


test_that("grouped LRT matches manual refits scored on the constant-free scale", {
  skip_if_not_installed("VGAM")

  dat <- make_wine_test_data()

  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  expect_identical(fit$stage1_engine, "vglm")
  expect_identical(fit$stage2_engine, "vglm")

  term <- "contact"
  f1_reduced <- stats::update(fit$stage1_formula, paste(". ~ . -", term))
  f2_reduced <- stats::update(fit$stage2_formula, paste(". ~ . -", term))

  fit1_reduced <- VGAM::vglm(
    f1_reduced,
    data = fit$data_stage1,
    family = VGAM::cumulative(link = "logitlink", parallel = TRUE, reverse = FALSE)
  )
  fit2_reduced <- VGAM::vglm(
    f2_reduced,
    data = fit$data_stage2,
    family = VGAM::multinomial(refLevel = 1)
  )

  # Grouped vglm log-likelihoods carry the multinomial constants, so they are
  # not comparable with `fit$logLik_stage*`; score the reduced fits the way the
  # package does, against the stored weighted counts.
  p1 <- .tsco_align_prob(
    VGAM::predictvglm(fit1_reduced, newdata = fit$data_stage1, type = "response"),
    .tsco_stage_levels(fit, 1L)
  )
  p2 <- .tsco_align_prob(
    VGAM::predictvglm(fit2_reduced, newdata = fit$data_stage2, type = "response"),
    .tsco_stage_levels(fit, 2L)
  )

  ll_reduced <- .tsco_count_loglik(fit$counts_stage1, p1) +
    .tsco_count_loglik(fit$counts_stage2, p2)

  lrt_manual <- 2 * (fit$logLik_stage1 + fit$logLik_stage2 - ll_reduced)

  # The native constants differ from the constant-free scale by a fixed
  # amount, so a naive native-scale statistic still agrees here (both fits
  # share the counts) -- but the *level* of the log-likelihood does not.
  expect_false(isTRUE(all.equal(
    as.numeric(fit1_reduced@criterion["loglikelihood"]),
    .tsco_count_loglik(fit$counts_stage1, p1),
    tolerance = 1e-3
  )))

  joint_tests <- summary(fit, joint_test = "LRT")$joint_tests

  expect_equal(joint_tests[term, "Chisq"], lrt_manual, tolerance = 1e-8)

  # And the same statistic as the individual-level representation.
  fit_ind <- tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )
  joint_ind <- summary(fit_ind, joint_test = "LRT")$joint_tests

  expect_equal(joint_tests[term, "Chisq"], joint_ind[term, "Chisq"], tolerance = 1e-6)
})
