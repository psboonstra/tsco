test_that("coef and vcov apply sign and naming transformations consistently", {
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

    # ------------------------------------------------------------------
    # Expected coefficients: sign correction plus reported names
    # ------------------------------------------------------------------

    b1_raw <- stats::coef(fit$fit_stage1)
    b2_raw <- stats::coef(fit$fit_stage2)

    s1 <- .tsco_sign_vector(
      fit$fit_stage1,
      kind = fit$stage1,
      engine = fit$stage1_engine
    )

    s2 <- .tsco_sign_vector(
      fit$fit_stage2,
      kind = fit$stage2,
      engine = fit$stage2_engine
    )

    b1_expected <- b1_raw * s1
    b2_expected <- b2_raw * s2

    map1 <- .tsco_stage_reported_names(
      fit$fit_stage1,
      kind = fit$stage1,
      engine = fit$stage1_engine,
      response_levels = fit$lower_levels,
      formula = fit$stage1_formula,
      data = fit$data_stage1,
      collapsed_label = fit$lower_collapsed_label,
      cutoff_level = fit$cutoff_level
    )

    map2 <- .tsco_stage_reported_names(
      fit$fit_stage2,
      kind = fit$stage2,
      engine = fit$stage2_engine,
      response_levels = fit$collapsed_levels,
      formula = fit$stage2_formula,
      data = fit$data_stage2,
      collapsed_label = fit$lower_collapsed_label,
      cutoff_level = fit$cutoff_level
    )

    names(b1_expected) <- .tsco_apply_reported_names(names(b1_raw), map1)
    names(b2_expected) <- .tsco_apply_reported_names(names(b2_raw), map2)

    expect_equal(
      coef(fit, stage = "stage1"),
      b1_expected
    )

    expect_equal(
      coef(fit, stage = "stage2"),
      b2_expected
    )

    # ------------------------------------------------------------------
    # Expected covariance: D V D plus reported names
    # ------------------------------------------------------------------

    V1_raw <- .tsco_stage_vcov(fit$fit_stage1)
    V2_raw <- .tsco_stage_vcov(fit$fit_stage2)

    s1_v <- .tsco_align_sign_vector(s1, V1_raw)
    s2_v <- .tsco_align_sign_vector(s2, V2_raw)

    V1_expected <- V1_raw * outer(s1_v, s1_v)
    V2_expected <- V2_raw * outer(s2_v, s2_v)

    rn1 <- .tsco_apply_reported_names(rownames(V1_expected), map1)
    rn2 <- .tsco_apply_reported_names(rownames(V2_expected), map2)

    dimnames(V1_expected) <- list(rn1, rn1)
    dimnames(V2_expected) <- list(rn2, rn2)

    expect_equal(
      vcov(fit, stage = "stage1"),
      V1_expected,
      tolerance = 1e-10
    )

    expect_equal(
      vcov(fit, stage = "stage2"),
      V2_expected,
      tolerance = 1e-10
    )
  }
})
