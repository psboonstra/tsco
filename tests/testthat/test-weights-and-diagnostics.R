# Frequency weights: sample-size and BIC convention, zero-weight rows, and the
# slope-only separation heuristic.

fit_wine <- function(data, rhs = "temp + contact", cutoff = NULL, ...) {
  dat <- make_wine_test_data()

  tsco(
    stats::as.formula(paste0("rating ~ ", rhs)),
    data = data,
    levels = dat$levels,
    cutoff_level = if (is.null(cutoff)) dat$cutoff_level else cutoff,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE,
    ...
  )
}


# --- weighted sample size and BIC -------------------------------------------

test_that("integer weights and row replication give the same n_weighted, nobs and BIC", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  set.seed(20260917)
  w <- sample(1:3, nrow(wine), replace = TRUE)
  wine_rep <- wine[rep(seq_len(nrow(wine)), times = w), , drop = FALSE]

  fit_w <- fit_wine(wine, weights = w)
  fit_r <- fit_wine(wine_rep)

  expect_true(fit_w$weighted)
  expect_false(fit_r$weighted)

  # Unweighted counts describe the rows supplied; weighted counts describe the
  # replicated sample the likelihood is on.
  expect_equal(fit_w$n_obs, nrow(wine))
  expect_equal(fit_w$n_weighted, sum(w))
  expect_equal(fit_r$n_obs, sum(w))
  expect_equal(fit_r$n_weighted, sum(w))

  lower <- wine$rating %in% dat$levels[1:3]
  expect_equal(fit_w$n_stage1_weighted, sum(w[lower]))
  expect_equal(fit_w$n_stage1_weighted, sum(fit_w$totals_stage1))

  expect_equal(attr(logLik(fit_w), "nobs"), sum(w))
  expect_equal(as.numeric(logLik(fit_w)), as.numeric(logLik(fit_r)), tolerance = 1e-6)
  expect_equal(BIC(fit_w), BIC(fit_r), tolerance = 1e-6)

  sw <- summary(fit_w, joint_test = "none")
  sr <- summary(fit_r, joint_test = "none")

  expect_equal(sw$component_fit$BIC, sr$component_fit$BIC, tolerance = 1e-6)
  expect_equal(sw$total_fit$BIC, sr$total_fit$BIC, tolerance = 1e-6)
  expect_equal(sw$component_fit$n_weighted, sr$component_fit$n_weighted)

  expect_output(print(sw), "Sum of weights")
  expect_output(print(sw), "frequency-weighted")
  expect_no_match(paste(utils::capture.output(print(sr)), collapse = "\n"), "Sum of weights")
})


test_that("grouped weights count the weighted total of each row's counts", {
  dat <- make_wine_test_data()
  w <- c(2, 1, 3, 1)

  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    weights = w,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  counts <- as.matrix(dat$grouped[dat$levels])
  expect_equal(fit$n_obs, sum(counts))
  expect_equal(fit$n_weighted, sum(w * rowSums(counts)))
  expect_equal(attr(logLik(fit), "nobs"), sum(w * rowSums(counts)))
})


test_that("without weights n_weighted equals n_obs and nothing about BIC changes", {
  dat <- make_wine_test_data()
  fit <- fit_wine(dat$individual)

  expect_equal(fit$n_weighted, fit$n_obs)
  expect_equal(fit$n_stage1_weighted, fit$n_stage1_obs)
  expect_equal(fit$n_zero_weight, 0L)
  expect_equal(
    BIC(fit),
    -2 * as.numeric(logLik(fit)) + log(fit$n_obs) * attr(logLik(fit), "df")
  )
})


# --- zero-weight rows -------------------------------------------------------

test_that("zero-weight rows are removed before fitting", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  set.seed(1)
  w <- sample(c(0, 1, 2), nrow(wine), replace = TRUE, prob = c(0.25, 0.5, 0.25))
  stopifnot(any(w == 0), any(w == 2))

  fit_w <- fit_wine(wine, weights = w, po_engine = "vglm")
  fit_sub <- fit_wine(wine[w > 0, , drop = FALSE], weights = w[w > 0], po_engine = "vglm")

  expect_equal(fit_w$n_zero_weight, sum(w == 0))
  expect_equal(fit_w$n_obs, sum(w > 0))
  expect_equal(fit_w$n_weighted, sum(w))
  expect_equal(coef(fit_w), coef(fit_sub), tolerance = 1e-6)
  expect_equal(fit_w$logLik, fit_sub$logLik, tolerance = 1e-6)
  expect_equal(fit_w$totals_stage1, fit_sub$totals_stage1)

  # predict() without newdata still covers every input row.
  p <- predict(fit_w)
  expect_equal(nrow(p), nrow(wine))
  expect_true(all(is.finite(p)))
  expect_equal(unname(rowSums(p)), rep(1, nrow(wine)), tolerance = 1e-8)
})


test_that("a factor level seen only in zero-weight rows establishes no support", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  wine$grp <- factor(
    ifelse(seq_len(nrow(wine)) <= 6, "U", ifelse(wine$temp == "cold", "A", "B")),
    levels = c("A", "B", "U")
  )
  w <- ifelse(wine$grp == "U", 0, 1)
  stopifnot(any(wine$rating[wine$grp == "U"] %in% dat$levels[1:3]))

  fit <- fit_wine(wine, rhs = "grp", weights = w, po_engine = "vglm")

  expect_equal(fit$n_zero_weight, 6L)
  expect_identical(fit$predictor_observed_levels$grp, c("A", "B"))
  expect_identical(fit$stage1_factor_levels$grp, c("A", "B"))
  expect_false(any(grepl("grpU", names(coef(fit)))))

  # Same fit as if the zero-weight rows had never been supplied.
  fit_ref <- fit_wine(wine[w > 0, , drop = FALSE], rhs = "grp", po_engine = "vglm")
  expect_equal(coef(fit), coef(fit_ref), tolerance = 1e-6)

  # The level cannot be predicted from newdata ...
  expect_error(
    predict(fit, newdata = data.frame(grp = factor("U", levels = c("A", "B", "U")))),
    "not observed in the fitting data"
  )

  # ... and predict(fit) returns those rows as all-NA, with a warning.
  expect_warning(p <- predict(fit), "zero weight")
  expect_equal(nrow(p), nrow(wine))
  expect_true(all(is.na(p[1:6, ])))
  expect_true(all(is.finite(p[-(1:6), ])))
  expect_identical(unname(attr(p, "unfitted")), w == 0)

  expect_warning(cls <- predict(fit, type = "class"), "zero weight")
  expect_true(all(is.na(cls[1:6])))
  expect_false(anyNA(cls[-(1:6)]))
})


test_that("an outcome level seen only in zero-weight rows is treated as unobserved", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  w <- ifelse(wine$rating == "2", 0, 1)

  expect_warning(
    fit <- fit_wine(wine, weights = w, po_engine = "vglm"),
    "declared in `levels` but have no observations"
  )

  expect_identical(fit$unobserved_levels, "2")
  expect_identical(names(fit$totals_stage1), c("1", "3"))

  p <- predict(fit)
  expect_equal(nrow(p), nrow(wine))
  expect_equal(unname(p[, "2"]), rep(0, nrow(wine)))
})


test_that("a grouped row with zero weight is removed", {
  dat <- make_wine_test_data()

  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    weights = c(0, 1, 1, 1),
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  fit_ref <- tsco(
    dat$grouped_formula,
    data = dat$grouped[-1, , drop = FALSE],
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  counts <- as.matrix(dat$grouped[dat$levels])
  expect_equal(fit$n_zero_weight, 1L)
  expect_equal(fit$n_groups, 3L)
  expect_equal(fit$n_obs, sum(counts[-1, ]))
  expect_equal(coef(fit), coef(fit_ref), tolerance = 1e-6)
  expect_equal(nrow(predict(fit)), nrow(dat$grouped))
})


# --- separation heuristic ---------------------------------------------------

test_that("the separation heuristic looks at slope rows only", {
  make_table <- function(est, se, names) {
    out <- cbind(Estimate = est, `Std. Error` = se, `z value` = est / se,
                 `Pr(>|z|)` = NA_real_)
    rownames(out) <- names
    out
  }

  fake_fit <- structure(list(fail = FALSE), class = "orm")

  # Unavailable threshold covariance is not separation.
  tab <- make_table(
    est = c(-1, 0.5, 2, 0.3),
    se = c(NA, 0.2, NA, 0.1),
    names = c("Y<=1", "Y<=2", "Y<=3", "x")
  )
  d <- .tsco_stage_diagnostics(fake_fit, engine = "orm", coef_table = tab)
  expect_false(d$separated)
  expect_length(d$messages, 0L)

  # A finite slope with a nonfinite standard error is.
  tab <- make_table(
    est = c(-1, 0.5, 8, 0.3),
    se = c(0.3, 0.2, NA, 0.1),
    names = c("Y<=1", "Y<=2", "x", "z")
  )
  d <- .tsco_stage_diagnostics(fake_fit, engine = "orm", coef_table = tab)
  expect_true(d$separated)
  expect_match(d$messages, "coefficient\\(s\\) x ")
  expect_no_match(d$messages, "Y<=")

  # So is the classical large-estimate, large-SE signature, but not on a
  # threshold, however extreme.
  tab <- make_table(
    est = c(-40, 0.5, 25, 0.3),
    se = c(300, 0.2, 400, 0.1),
    names = c("Y<=1", "Y<=2", "x", "z")
  )
  d <- .tsco_stage_diagnostics(fake_fit, engine = "orm", coef_table = tab)
  expect_true(d$separated)
  expect_match(d$messages, "coefficient\\(s\\) x ")

  # Multinomial intercepts are thresholds too.
  tab <- make_table(
    est = c(-40, 12, 0.3),
    se = c(NA, 90, 0.1),
    names = c("(Intercept):2", "x:2", "x:3")
  )
  d <- .tsco_stage_diagnostics(fake_fit, engine = "vglm", coef_table = tab)
  expect_true(d$separated)
  expect_match(d$messages, "coefficient\\(s\\) x:2 ")
})


test_that(".tsco_is_threshold_name recognises reported threshold names", {
  expect_identical(
    .tsco_is_threshold_name(c("Y<=1", "Y<4", "(Intercept):2", "x", "x:2",
                              "stage1:Y<=2", "stage2:tempwarm", "log(z)")),
    c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)
  )
})


# --- binary multinomial names -----------------------------------------------

test_that("a two-category multinomial stage reports category-labelled names", {
  dat <- make_wine_test_data()

  # Cutoff at the last level makes the stage-2 collapsed outcome binary.
  fit <- tsco(
    rating ~ temp,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = "5",
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  expect_identical(.tsco_stage_levels(fit, 2L), c(fit$lower_collapsed_label, "5"))

  s2 <- grep("^stage2:", names(coef(fit)), value = TRUE)
  expect_identical(s2, c("stage2:(Intercept):5", "stage2:tempwarm:5"))
  expect_identical(rownames(vcov(fit, stage = "stage2")), c("(Intercept):5", "tempwarm:5"))
  expect_identical(
    rownames(summary(fit, joint_test = "none")$coef_stage2),
    c("(Intercept):5", "tempwarm:5")
  )
})


# --- no-intercept formulas ---------------------------------------------------

test_that("a formula without an intercept is an error that points at relevel()", {
  dat <- make_wine_test_data()

  for (rhs in c("temp - 1", "0 + temp", "temp + contact - 1")) {
    expect_error(
      fit_wine(dat$individual, rhs = rhs),
      "must include an intercept.*relevel"
    )
  }

  expect_error(
    tsco(
      update(dat$grouped_formula, . ~ . - 1),
      data = dat$grouped,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = "po",
      stage2 = "multinomial",
      warn_degenerate = FALSE
    ),
    "must include an intercept"
  )

  # The recommended route works and changes only the reference level.
  wine <- dat$individual
  wine$temp <- stats::relevel(factor(wine$temp, ordered = FALSE), ref = "warm")
  fit <- fit_wine(wine, rhs = "temp")
  expect_true("stage1:tempcold" %in% names(coef(fit)))
})


# --- uniform completion and unobserved lower levels ---------------------------

test_that("uniform completion puts no mass on a declared but unobserved lower level", {
  # Sparse stage-1 fixture: level "U" of grp never appears below the cutoff,
  # and lower level "b" is declared but never observed at all.
  grouped <- data.frame(
    grp = factor(c("A", "B", "U"), levels = c("A", "B", "U")),
    a = c(18, 8, 0),
    b = c(0, 0, 0),
    c = c(7, 17, 0),
    d = c(8, 6, 18),
    e = c(4, 9, 12)
  )

  expect_warning(
    fit <- tsco(
      cbind(a, b, c, d, e) ~ grp,
      data = grouped,
      levels = c("a", "b", "c", "d", "e"),
      cutoff_level = "d",
      stage1 = "po",
      stage2 = "po",
      po_engine = "vglm",
      warn_degenerate = FALSE
    ),
    "declared in `levels` but have no observations"
  )

  expect_identical(fit$unobserved_levels, "b")
  expect_identical(.tsco_stage_levels(fit, 1L), c("a", "c"))

  nd <- data.frame(grp = factor("U", levels = c("A", "B", "U")))

  for (method in c("uniform", "empirical")) {
    expect_warning(
      p <- predict(fit, newdata = nd, unsupported_stage1 = method),
      "not represented in the stage-1 fitting data"
    )
    expect_equal(unname(p[, "b"]), 0)
    expect_equal(unname(rowSums(p)), 1, tolerance = 1e-10)
  }

  p_u <- suppressWarnings(predict(fit, newdata = nd, unsupported_stage1 = "uniform"))
  lm <- unname(attr(p_u, "lower_mass"))
  expect_equal(unname(p_u[, c("a", "b", "c")]), c(lm / 2, 0, lm / 2), tolerance = 1e-10)

  p_e <- suppressWarnings(predict(fit, newdata = nd, unsupported_stage1 = "empirical"))
  emp <- c(a = 26, c = 24) / 50
  expect_equal(unname(p_e[, c("a", "c")]), unname(lm * emp), tolerance = 1e-10)
})


test_that("uniform completion is uniform over the observed lower levels", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  # grp "U" appears only above the cutoff, so it is unsupported in stage 1;
  # rating 2 is declared but never observed.
  wine <- wine[wine$rating != "2", , drop = FALSE]
  wine$grp <- factor(ifelse(wine$temp == "cold", "A", "B"), levels = c("A", "B", "U"))
  wine$grp[wine$rating %in% c("4", "5")][1:4] <- "U"
  stopifnot(!any(wine$grp == "U" & wine$rating %in% c("1", "3")))

  fit <- suppressWarnings(fit_wine(wine, rhs = "grp", po_engine = "vglm"))
  expect_identical(fit$unobserved_levels, "2")
  expect_identical(.tsco_stage_levels(fit, 1L), c("1", "3"))

  nd <- data.frame(grp = factor("U", levels = c("A", "B", "U")))
  p <- suppressWarnings(predict(fit, newdata = nd, unsupported_stage1 = "uniform"))
  lm <- unname(attr(p, "lower_mass"))

  expect_equal(unname(p[, c("1", "2", "3")]), c(lm / 2, 0, lm / 2), tolerance = 1e-10)
  expect_equal(unname(rowSums(p)), 1, tolerance = 1e-10)
})
