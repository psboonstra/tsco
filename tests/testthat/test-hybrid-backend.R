test_that("individual-level PO stages use rms::orm", {
  skip_if_not_installed("rms")
  skip_if_not_installed("VGAM")

  set.seed(6819)
  dat <- data.frame(
    x = stats::rnorm(400),
    y = ordered(
      sample(c("a", "b", "c", "d"), 400, replace = TRUE),
      levels = c("a", "b", "c", "d")
    )
  )

  fit_popo <- tsco(
    y ~ x,
    data = dat,
    cutoff_level = "c",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_identical(fit_popo$stage1_engine, "orm")
  expect_identical(fit_popo$stage2_engine, "orm")
  expect_s3_class(fit_popo$fit_stage1, "orm")
  expect_s3_class(fit_popo$fit_stage2, "orm")
  expect_false("po.reverse" %in% names(formals(tsco)))
  expect_identical(intersect("...", names(formals(tsco))), character())
  expect_false("po.reverse" %in% names(fit_popo))
  expect_false("vglm_args" %in% names(fit_popo))
  expect_identical(intersect("fit_args", names(fit_popo)), character())
  expect_identical(fit_popo$stage1_args, list())
  expect_identical(fit_popo$stage2_args, list())
  expect_equal(unname(rowSums(predict(fit_popo))), rep(1, nrow(dat)), tolerance = 1e-10)
  expect_true(all(is.finite(summary(fit_popo, joint_test = "LRT")$joint_tests$Chisq)))
  expect_identical(rownames(summary(fit_popo, joint_test = "Wald")$joint_tests), "x")

  fit_pomr <- tsco(
    y ~ x,
    data = dat,
    cutoff_level = "c",
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  expect_identical(fit_pomr$stage1_engine, "orm")
  expect_identical(fit_pomr$stage2_engine, "vglm")
  expect_s3_class(fit_pomr$fit_stage1, "orm")
  expect_s4_class(fit_pomr$fit_stage2, "vglm")

  fit_orm_popo <- tsco(
    y ~ x,
    data = dat,
    cutoff_level = "c",
    stage1 = "po",
    stage2 = "po",
    po_engine = "orm",
    warn_degenerate = FALSE
  )

  expect_identical(fit_orm_popo$po_engine, "orm")
  expect_identical(fit_orm_popo$stage1_engine, "orm")
  expect_identical(fit_orm_popo$stage2_engine, "orm")
  expect_s3_class(fit_orm_popo$fit_stage1, "orm")
  expect_s3_class(fit_orm_popo$fit_stage2, "orm")
  expect_equal(predict(fit_orm_popo), predict(fit_popo), tolerance = 1e-10)

  fit_vglm_popo <- tsco(
    y ~ x,
    data = dat,
    cutoff_level = "c",
    stage1 = "po",
    stage2 = "po",
    po_engine = "vglm",
    warn_degenerate = FALSE
  )

  expect_identical(fit_vglm_popo$po_engine, "vglm")
  expect_identical(fit_vglm_popo$stage1_engine, "vglm")
  expect_identical(fit_vglm_popo$stage2_engine, "vglm")
  expect_s4_class(fit_vglm_popo$fit_stage1, "vglm")
  expect_s4_class(fit_vglm_popo$fit_stage2, "vglm")

  orm_stage1 <- coef(fit_popo, stage = "stage1")
  vglm_stage1 <- coef(fit_vglm_popo, stage = "stage1")
  orm_stage1_thresholds <- orm_stage1[
    seq_len(fit_popo$fit_stage1$non.slopes)
  ]
  vglm_stage1_thresholds <- vglm_stage1[
    grepl("(Intercept)", names(vglm_stage1), fixed = TRUE)
  ]
  expect_equal(
    unname(orm_stage1_thresholds),
    unname(vglm_stage1_thresholds),
    tolerance = 1e-5
  )

  orm_stage2 <- coef(fit_popo, stage = "stage2")
  vglm_stage2 <- coef(fit_vglm_popo, stage = "stage2")
  orm_stage2_thresholds <- orm_stage2[
    seq_len(fit_popo$fit_stage2$non.slopes)
  ]
  vglm_stage2_thresholds <- vglm_stage2[
    grepl("(Intercept)", names(vglm_stage2), fixed = TRUE)
  ]
  expect_equal(
    unname(orm_stage2_thresholds),
    unname(vglm_stage2_thresholds),
    tolerance = 1e-5
  )

  expect_equal(predict(fit_vglm_popo), predict(fit_popo), tolerance = 1e-6)
  expect_true(all(is.finite(summary(fit_vglm_popo, joint_test = "LRT")$joint_tests$Chisq)))
})


test_that("many individual PO levels use rms::orm", {
  skip_if_not_installed("rms")

  levels <- c(paste0("L", seq_len(40)), "terminal")
  dat <- data.frame(
    x = rep(seq(-1, 1, length.out = 10), each = 41),
    y = ordered(rep(levels, times = 10), levels = levels)
  )

  fit <- tsco(
    y ~ x,
    data = dat,
    cutoff_level = "terminal",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_identical(fit$stage1_engine, "orm")
  expect_s3_class(fit$fit_stage1, "orm")
  expect_equal(fit$fit_stage1$non.slopes, 39L)

  stage1_coef <- coef(fit, stage = "stage1")
  V_full <- vcov(fit, stage = "stage1")
  expect_equal(dim(V_full), rep(length(stage1_coef), 2L))
  expect_identical(rownames(V_full), names(stage1_coef))
  expect_identical(colnames(V_full), names(stage1_coef))
  expect_equal(sum(!is.finite(V_full)), 0L)

  V_mid <- stats::vcov(fit$fit_stage1)
  reduced_orm <- fit$fit_stage1
  reduced_orm$var <- V_mid
  V_reduced <- .tsco_stage_vcov(reduced_orm)
  unavailable <- setdiff(names(stats::coef(reduced_orm)), rownames(V_mid))

  expect_equal(dim(V_reduced), rep(length(stage1_coef), 2L))
  expect_equal(
    V_reduced[rownames(V_mid), colnames(V_mid), drop = FALSE],
    as.matrix(V_mid)
  )
  expect_equal(
    unname(diag(V_reduced)[unavailable]),
    rep(NA_real_, length(unavailable))
  )
  expect_equal(sum(!is.finite(V_reduced["x", "x", drop = FALSE])), 0L)

  expect_equal(unname(rowSums(predict(fit))), rep(1, nrow(dat)), tolerance = 1e-10)
})


test_that("grouped PO stages retain the VGAM backend", {
  skip_if_not_installed("VGAM")

  dat <- data.frame(
    x = c(0, 1, 2),
    a = c(12, 8, 5),
    b = c(6, 8, 9),
    c = c(4, 6, 10),
    d = c(2, 4, 6)
  )

  fit <- tsco(
    cbind(a, b, c, d) ~ x,
    data = dat,
    cutoff_level = "c",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_identical(fit$stage1_engine, "vglm")
  expect_identical(fit$stage2_engine, "vglm")
  expect_s4_class(fit$fit_stage1, "vglm")
  expect_s4_class(fit$fit_stage2, "vglm")

  expect_snapshot(
    tsco(
      cbind(a, b, c, d) ~ x,
      data = dat,
      cutoff_level = "c",
      stage1 = "po",
      stage2 = "po",
      po_engine = "orm",
      warn_degenerate = FALSE
    ),
    error = TRUE
  )
})
