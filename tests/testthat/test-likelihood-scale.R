# Every log-likelihood entering a likelihood-ratio statistic must be on the
# package's constant-free scale. Grouped-count `VGAM::vglm()` fits carry the
# multinomial combinatorial constants in their own log-likelihood; `rms::orm()`
# and the closed-form intercept-only fit do not. The tests here pin the LRT to
# the individual-level answer, which is unaffected by the constants.

fit_wine_pair <- function(rhs, stage1, stage2, weights = NULL) {
  dat <- make_wine_test_data()

  grouped_formula <- stats::as.formula(
    paste0(.tsco_deparse(dat$grouped_formula[[2L]]), " ~ ", rhs)
  )
  individual_formula <- stats::as.formula(paste0("rating ~ ", rhs))

  list(
    grouped = tsco(
      grouped_formula,
      data = dat$grouped,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = stage1,
      stage2 = stage2,
      warn_degenerate = FALSE
    ),
    individual = tsco(
      individual_formula,
      data = dat$individual,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = stage1,
      stage2 = stage2,
      warn_degenerate = FALSE
    )
  )
}

expect_same_lrt <- function(fits) {
  lg <- summary(fits$grouped, joint_test = "LRT")$joint_tests
  li <- summary(fits$individual, joint_test = "LRT")$joint_tests

  # The grouped and individual fits take different IRLS paths, so agreement
  # is to fitter precision (about 1e-7 in the statistic), not machine
  # precision. A likelihood-constant mismatch would be off by tens.
  expect_equal(rownames(lg), rownames(li))
  expect_equal(lg$df, li$df)
  expect_equal(lg$Chisq, li$Chisq, tolerance = 1e-6)
  expect_equal(lg[["Pr(>Chisq)"]], li[["Pr(>Chisq)"]], tolerance = 1e-5)

  invisible(list(grouped = lg, individual = li))
}


test_that("a single-predictor grouped LRT agrees with the individual-level LRT", {
  # This is the manuscript's own `y ~ treatment` design. Before the fix the
  # grouped statistic was 162.2 against a correct 24.2.
  fits <- fit_wine_pair("temp", stage1 = "po", stage2 = "po")
  out <- expect_same_lrt(fits)

  expect_true(all(is.finite(out$grouped$Chisq)))
  expect_lt(out$grouped$Chisq, 50)
})


test_that("single-predictor grouped LRTs agree for multinomial stages too", {
  fits <- fit_wine_pair("temp", stage1 = "multinomial", stage2 = "multinomial")
  expect_same_lrt(fits)

  fits <- fit_wine_pair("temp", stage1 = "po", stage2 = "multinomial")
  expect_same_lrt(fits)
})


test_that("multi-term grouped LRTs still agree with the individual-level LRT", {
  fits <- fit_wine_pair("temp + contact", stage1 = "po", stage2 = "po")
  expect_same_lrt(fits)

  fits <- fit_wine_pair("temp + contact", stage1 = "po", stage2 = "multinomial")
  expect_same_lrt(fits)
})


test_that("the single-predictor LRT equals the closed-form statistic", {
  fits <- fit_wine_pair("temp", stage1 = "po", stage2 = "po")

  for (fit in fits) {
    null1 <- .tsco_null_stage_fit_stats(fit$totals_stage1)
    null2 <- .tsco_null_stage_fit_stats(fit$totals_stage2)

    manual <- 2 * (
      fit$logLik_stage1 + fit$logLik_stage2 - null1$loglik - null2$loglik
    )

    lrt <- summary(fit, joint_test = "LRT")$joint_tests

    expect_equal(lrt$Chisq, manual, tolerance = 1e-8)
  }
})


test_that("frequency weights and row replication give the same LRT", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  set.seed(20260917)
  w <- sample(1:3, nrow(wine), replace = TRUE)
  wine_rep <- wine[rep(seq_len(nrow(wine)), times = w), , drop = FALSE]

  for (rhs in c("temp", "temp + contact")) {
    f <- stats::as.formula(paste0("rating ~ ", rhs))

    fit_w <- tsco(
      f,
      data = wine,
      weights = w,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = "po",
      stage2 = "po",
      warn_degenerate = FALSE
    )

    fit_r <- tsco(
      f,
      data = wine_rep,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = "po",
      stage2 = "po",
      warn_degenerate = FALSE
    )

    lw <- summary(fit_w, joint_test = "LRT")$joint_tests
    lr <- summary(fit_r, joint_test = "LRT")$joint_tests

    expect_equal(lw$df, lr$df)
    expect_equal(lw$Chisq, lr$Chisq, tolerance = 1e-6)
  }
})


test_that("stored stage count matrices reproduce the stored log-likelihoods", {
  fits <- fit_wine_pair("temp + contact", stage1 = "po", stage2 = "multinomial")

  for (fit in fits) {
    expect_equal(colnames(fit$counts_stage1), .tsco_stage_levels(fit, 1L))
    expect_equal(colnames(fit$counts_stage2), .tsco_stage_levels(fit, 2L))
    expect_equal(nrow(fit$counts_stage1), nrow(fit$data_stage1))
    expect_equal(nrow(fit$counts_stage2), nrow(fit$data_stage2))
    expect_equal(colSums(fit$counts_stage1), fit$totals_stage1)
    expect_equal(colSums(fit$counts_stage2), fit$totals_stage2)

    p1 <- .tsco_align_prob(
      .tsco_predict_stage_prob(
        fit$fit_stage1, fit$stage1_engine, fit$data_stage1,
        levels = .tsco_stage_levels(fit, 1L)
      ),
      .tsco_stage_levels(fit, 1L)
    )

    expect_equal(.tsco_count_loglik(fit$counts_stage1, p1), fit$logLik_stage1)
  }
})


# --- empirical completion of unsupported stage-1 rows -----------------------

test_that("the empirical lower-category distribution gives an unobserved level zero mass", {
  dat <- make_wine_test_data()

  # Rating 2 is a declared lower level with no observations.
  grouped <- dat$grouped
  grouped[["2"]] <- 0L
  individual <- dat$individual[dat$individual$rating != "2", , drop = FALSE]

  fit_g <- suppressWarnings(tsco(
    dat$grouped_formula,
    data = grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  ))

  fit_i <- suppressWarnings(tsco(
    dat$individual_formula,
    data = individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  ))

  expected <- table(factor(individual$rating[individual$rating %in% c("1", "2", "3")],
                           levels = c("1", "2", "3")))
  expected <- as.numeric(expected) / sum(expected)
  names(expected) <- c("1", "2", "3")

  expect_no_warning(emp_g <- .tsco_empirical_lower_probs(fit_g))
  expect_no_warning(emp_i <- .tsco_empirical_lower_probs(fit_i))

  expect_equal(emp_g, expected)
  expect_equal(emp_i, expected)
  expect_identical(unname(emp_g["2"]), 0)
})


test_that("the empirical lower-category distribution honours case weights", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  set.seed(1)
  w <- sample(1:3, nrow(wine), replace = TRUE)

  fit_w <- tsco(
    rating ~ temp + contact,
    data = wine,
    weights = w,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  lower <- wine$rating %in% dat$levels[1:3]
  expected <- tapply(w[lower], factor(wine$rating[lower], levels = dat$levels[1:3]), sum)
  expected <- as.numeric(expected) / sum(expected)
  names(expected) <- dat$levels[1:3]

  expect_equal(.tsco_empirical_lower_probs(fit_w), expected)
})


# --- basis terms in the formula --------------------------------------------

test_that("basis terms such as poly() predict identically from raw newdata", {
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$z <- seq_len(nrow(wine))

  fit <- tsco(
    rating ~ poly(z, 2) + temp,
    data = wine,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  # The stored terms must carry the fitted basis constants; otherwise
  # `poly(z, 2)` is rebuilt from `newdata` alone and the predictions are
  # silently wrong (and a two-row `newdata` is an error).
  expect_false(is.null(attr(fit$terms_rhs, "predvars")))

  rows <- c(10L, 20L, 30L)
  p_fit <- predict(fit)
  p_new <- predict(fit, newdata = wine[rows, c("z", "temp")])

  expect_equal(unclass(p_new)[, ], p_fit[rows, ], tolerance = 1e-10,
               ignore_attr = TRUE)

  p_two <- predict(fit, newdata = wine[1:2, c("z", "temp")])
  expect_equal(nrow(p_two), 2L)
  expect_equal(unname(rowSums(p_two)), c(1, 1), tolerance = 1e-8)
})


test_that("scale() in the formula predicts identically from raw newdata", {
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$z <- seq_len(nrow(wine))

  fit <- tsco(
    rating ~ scale(z) + temp,
    data = wine,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  rows <- c(5L, 40L)
  p_fit <- predict(fit)
  p_new <- predict(fit, newdata = wine[rows, c("z", "temp")])

  expect_equal(unclass(p_new)[, ], p_fit[rows, ], tolerance = 1e-10,
               ignore_attr = TRUE)
})
