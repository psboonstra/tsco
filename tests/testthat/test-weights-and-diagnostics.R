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


# --- weight alignment after missing-data removal -----------------------------

test_that("weights align with the model frame under custom and misleading row names", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  set.seed(7)
  w <- sample(1:3, nrow(wine), replace = TRUE)
  wine$contact[c(3, 40)] <- NA
  complete <- !is.na(wine$contact)

  reference <- fit_wine(wine[complete, , drop = FALSE], weights = w[complete], po_engine = "vglm")

  # Character row names used to be a hard error.
  wine_chr <- wine
  rownames(wine_chr) <- paste0("patient_", seq_len(nrow(wine)))
  fit_chr <- fit_wine(wine_chr, weights = w, po_engine = "vglm")
  expect_equal(coef(fit_chr), coef(reference), tolerance = 1e-8)
  expect_equal(fit_chr$n_weighted, sum(w[complete]))

  # Numeric-looking row names that are not row positions used to select the
  # wrong weights without any error.
  wine_num <- wine
  rownames(wine_num) <- as.character(seq_len(nrow(wine)) + 1000L)
  fit_num <- fit_wine(wine_num, weights = w, po_engine = "vglm")
  expect_equal(coef(fit_num), coef(reference), tolerance = 1e-8)
  expect_equal(fit_num$n_weighted, sum(w[complete]))

  # One weight per complete case is still accepted.
  fit_cc <- fit_wine(wine, weights = w[complete], po_engine = "vglm")
  expect_equal(coef(fit_cc), coef(reference), tolerance = 1e-8)

  # Any other length is an error.
  expect_error(fit_wine(wine, weights = w[-1]), "one entry per row of `data`")
})


test_that("weights that are positive only on rows removed by na.action are an error", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  w <- rep(0, nrow(wine))
  w[1:3] <- 1
  wine$contact[1:3] <- NA

  expect_error(
    fit_wine(wine, weights = w),
    "No positive-weight observations remain"
  )
})


# --- offsets are rejected ------------------------------------------------------

test_that("an offset in the formula is a clear error", {
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$os <- seq_len(nrow(wine)) / 100

  expect_error(fit_wine(wine, rhs = "temp + offset(os)"), "Offsets are not supported")

  grouped <- dat$grouped
  grouped$n_total <- rowSums(grouped[dat$levels])

  expect_error(
    tsco(
      update(dat$grouped_formula, . ~ . + offset(log(n_total))),
      data = grouped,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = "po",
      stage2 = "po",
      warn_degenerate = FALSE
    ),
    "Offsets are not supported"
  )

  expect_length(.tsco_stage_rhs_labels(stats::terms(y ~ a + b:c)), 2L)
})


# --- coefficient-table fallback keeps reported names --------------------------

test_that("the estimate-only coefficient table uses reported names", {
  dat <- make_wine_test_data()
  fit <- fit_wine(dat$individual)

  reported <- .tsco_stage_reported_names(
    fit$fit_stage1,
    kind = fit$stage1,
    engine = fit$stage1_engine,
    response_levels = .tsco_stage_levels(fit, 1L),
    formula = fit$stage1_formula,
    data = fit$data_stage1,
    collapsed_label = fit$lower_collapsed_label,
    cutoff_level = fit$cutoff_level
  )

  full <- .tsco_stage_coef_table(
    fit$fit_stage1, kind = fit$stage1, engine = fit$stage1_engine,
    reported_names = reported
  )

  # Simulate a fit whose covariance cannot be extracted.
  local_mocked_bindings(.tsco_stage_vcov = function(fit) NULL, .package = "tsco")

  fallback <- .tsco_stage_coef_table(
    fit$fit_stage1, kind = fit$stage1, engine = fit$stage1_engine,
    reported_names = reported
  )

  expect_identical(colnames(fallback), "Estimate")
  expect_identical(rownames(fallback), rownames(full))
  expect_equal(fallback[, "Estimate"], full[, "Estimate"])
  expect_true(any(grepl("^Y<=", rownames(fallback))))
})


# --- rms's weights note is not passed on ------------------------------------

test_that("weighted orm fits do not pass on rms's validate()/bootcov() weights note", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  set.seed(11)
  w <- sample(1:2, nrow(wine), replace = TRUE)

  seen <- character()
  fit <- withCallingHandlers(
    {
      f <- fit_wine(wine, weights = w)
      summary(f, joint_test = "LRT")
      f
    },
    warning = function(cnd) {
      seen <<- c(seen, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(any(grepl("weights are ignored in model validation", seen, fixed = TRUE)))
  expect_true(fit$weighted)
})


# --- nobs(), required data, missing predictors in newdata ---------------------

test_that("nobs() returns the frequency-weighted sample size", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  fit <- fit_wine(wine)
  expect_equal(nobs(fit), nrow(wine))
  expect_equal(nobs(fit), attr(logLik(fit), "nobs"))

  set.seed(5)
  w <- sample(1:3, nrow(wine), replace = TRUE)
  fit_w <- fit_wine(wine, weights = w)
  expect_equal(nobs(fit_w), sum(w))
  expect_equal(nobs(fit_w), attr(logLik(fit_w), "nobs"))
  expect_equal(BIC(fit_w), -2 * as.numeric(logLik(fit_w)) + log(nobs(fit_w)) * attr(logLik(fit_w), "df"))
})


test_that("`data` is required and must be a data frame", {
  dat <- make_wine_test_data()
  rating <- dat$individual$rating
  temp <- dat$individual$temp

  expect_error(
    tsco(rating ~ temp, levels = dat$levels, cutoff_level = dat$cutoff_level,
         stage1 = "po", stage2 = "po", warn_degenerate = FALSE),
    "`data` must be supplied"
  )

  expect_error(
    tsco(rating ~ temp, data = as.list(dat$individual), levels = dat$levels,
         cutoff_level = dat$cutoff_level, stage1 = "po", stage2 = "po",
         warn_degenerate = FALSE),
    "must be a data frame"
  )
})


test_that("a missing predictor value in newdata gives an NA row, not a failure", {
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$z <- as.numeric(seq_len(nrow(wine)))

  nd <- data.frame(
    z = c(10, NA, 30, 40),
    temp = factor(c("cold", "warm", NA, "warm"), levels = levels(wine$temp))
  )
  nd_ok <- nd[c(1L, 4L), , drop = FALSE]

  for (engine in c("orm", "vglm")) {
    for (stage2 in c("po", "multinomial")) {
      fit <- tsco(
        rating ~ z + temp,
        data = wine,
        levels = dat$levels,
        cutoff_level = dat$cutoff_level,
        stage1 = "po",
        stage2 = stage2,
        po_engine = engine,
        warn_degenerate = FALSE
      )

      expect_warning(p <- predict(fit, newdata = nd), "missing or non-finite")

      expect_equal(nrow(p), 4L)
      expect_true(all(is.na(p[2:3, ])))
      expect_true(all(is.finite(p[c(1L, 4L), ])))
      expect_identical(unname(attr(p, "incomplete")), c(FALSE, TRUE, TRUE, FALSE))

      # Complete rows are unaffected by the presence of incomplete ones.
      p_ok <- predict(fit, newdata = nd_ok)
      expect_equal(unclass(p[c(1L, 4L), ]), unclass(p_ok)[, ], ignore_attr = TRUE)

      # An NA in a factor is a missing value, not an unsupported level.
      expect_false(any(attr(p, "unsupported_stage1")))

      expect_warning(cls <- predict(fit, newdata = nd, type = "class"), "missing or non-finite")
      expect_true(all(is.na(cls[2:3])))
      expect_false(anyNA(cls[c(1L, 4L)]))
    }
  }
})


# --- weights travel through model.frame() -------------------------------------

test_that("weights stay aligned under a custom na.action that records row names", {
  dat <- make_wine_test_data()
  wine <- dat$individual

  # A deliberately non-standard na.action: it removes incomplete rows but
  # records them by row *name*, not position, so generic `naresid()`-style
  # restoration would misread it. The package must not depend on how an
  # na.action records omissions, only on what it returns. The data's row
  # names are chosen not to be row positions.
  na_by_name <- function(object, ...) {
    keep <- stats::complete.cases(object)
    out <- object[keep, , drop = FALSE]
    omit <- as.integer(rownames(object)[!keep])
    names(omit) <- rownames(object)[!keep]
    class(omit) <- "omit"
    attr(out, "na.action") <- omit
    out
  }

  rownames(wine) <- as.character(seq_len(nrow(wine)) + 1000L)
  set.seed(9)
  w <- sample(1:3, nrow(wine), replace = TRUE)
  wine$contact[c(5, 30)] <- NA
  complete <- !is.na(wine$contact)

  reference <- fit_wine(wine[complete, , drop = FALSE], weights = w[complete], po_engine = "vglm")

  fit <- fit_wine(wine, weights = w, na.action = na_by_name, po_engine = "vglm")
  expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
  expect_equal(fit$n_weighted, sum(w[complete]))
  expect_equal(fit$n_obs, sum(complete))

  # The model frame's (weights) column must not leak into the predictors.
  expect_false("(weights)" %in% names(fit$predict_data))
  expect_false("(weights)" %in% names(fit$data_stage2))
  expect_false(any(grepl("weights", names(coef(fit)))))

  # na.exclude behaves like na.omit here.
  fit_ex <- fit_wine(wine, weights = w, na.action = stats::na.exclude, po_engine = "vglm")
  expect_equal(coef(fit_ex), coef(reference), tolerance = 1e-8)
})


test_that("the weight-length error names both accepted lengths", {
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$contact[1] <- NA

  expect_error(fit_wine(wine, weights = rep(1, 10)), "one entry per row of `data` or one per complete case")
})


# --- incomplete or non-finite fitting data are an early, clear error ----------

test_that("NA or non-finite predictors surviving na.action are a clear error", {
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$x <- as.numeric(seq_len(nrow(wine)))

  # na.pass leaves an NA factor in the frame.
  w_na <- wine
  w_na$temp[5] <- NA
  expect_error(
    fit_wine(w_na, rhs = "x + temp", na.action = stats::na.pass),
    "`temp`.*missing or non-finite"
  )

  # log(0) is -Inf after the formula is evaluated; na.omit does not remove it.
  w_inf <- wine
  w_inf$x[3] <- 0
  expect_error(
    fit_wine(w_inf, rhs = "log(x) + temp"),
    "`log\\(x\\)`.*missing or non-finite"
  )

  for (engine in c("orm", "vglm")) {
    expect_error(
      fit_wine(w_inf, rhs = "log(x) + temp", po_engine = engine),
      "missing or non-finite"
    )
  }

  # A missing response under na.pass is also caught.
  w_y <- wine
  w_y$rating[7] <- NA
  expect_error(
    fit_wine(w_y, rhs = "x + temp", na.action = stats::na.pass),
    "response contains missing values"
  )

  # na.fail is honoured as-is.
  expect_error(fit_wine(w_na, rhs = "x + temp", na.action = stats::na.fail), "missing values")
})


# --- non-finite transformed predictors in newdata ------------------------------

test_that("non-finite values after transformation in newdata give NA rows", {
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$x <- as.numeric(seq_len(nrow(wine)))

  nd <- data.frame(
    x = c(1, 0, NA, Inf, 20),
    temp = factor("cold", levels = levels(wine$temp))
  )

  for (engine in c("orm", "vglm")) {
    fit <- fit_wine(wine, rhs = "log(x) + temp", po_engine = engine)

    expect_warning(p <- predict(fit, newdata = nd), "missing or non-finite")
    expect_identical(unname(attr(p, "incomplete")), c(FALSE, TRUE, TRUE, TRUE, FALSE))
    expect_true(all(is.na(p[2:4, ])))
    expect_true(all(is.finite(p[c(1L, 5L), ])))

    p_ok <- predict(fit, newdata = nd[c(1L, 5L), , drop = FALSE])
    expect_equal(unclass(p[c(1L, 5L), ]), unclass(p_ok)[, ], ignore_attr = TRUE)

    st <- suppressWarnings(predict(fit, newdata = nd, type = "stage"))
    expect_identical(st$incomplete, attr(p, "incomplete"))
    expect_identical(st$unfitted, rep(FALSE, 5L))
    expect_true(all(is.na(st$prob[2:4, ])))
  }

  # Matrix-valued basis columns are screened cell-wise too.
  fit_p <- fit_wine(wine, rhs = "poly(x, 2) + temp")
  nd_p <- data.frame(x = c(5, NA, 10), temp = factor("warm", levels = levels(wine$temp)))
  expect_warning(pp <- predict(fit_p, newdata = nd_p), "missing or non-finite")
  expect_identical(unname(attr(pp, "incomplete")), c(FALSE, TRUE, FALSE))
  expect_true(all(is.finite(pp[c(1L, 3L), ])))
})


# --- weights do not leak into the stored terms ----------------------------------

test_that("passing weights through model.frame() leaves terms_rhs formula-only", {
  # `model.frame(weights = )` adds a `(weights)` column but must not add a
  # weights variable to the terms object; if it ever did, raw-`newdata`
  # prediction (which has no weights) would fail or, worse, evaluate wrongly.
  dat <- make_wine_test_data()
  wine <- dat$individual
  wine$x <- as.numeric(seq_len(nrow(wine)))

  set.seed(2)
  w <- sample(1:3, nrow(wine), replace = TRUE)
  wine$temp[c(2, 9)] <- NA

  for (rhs in c("poly(x, 2) + temp", "log(x) + temp", "scale(x) + contact")) {
    fit <- fit_wine(wine, rhs = rhs, weights = w, po_engine = "vglm")
    tt <- fit$terms_rhs

    vars <- vapply(as.list(attr(tt, "variables"))[-1L], deparse, character(1L))
    pv <- vapply(as.list(attr(tt, "predvars"))[-1L], function(v) deparse(v)[1L], character(1L))

    expect_false(any(grepl("weights", vars, fixed = TRUE)))
    expect_false(any(grepl("weights", attr(tt, "term.labels"), fixed = TRUE)))
    expect_false(any(grepl("weights", pv, fixed = TRUE)))
    expect_length(vars, 2L)
    expect_false("(weights)" %in% names(fit$predict_data))

    # Raw newdata (no weights column) predicts exactly what predict(fit) does.
    # predict(fit) is indexed by model-frame row, i.e. after na.omit, which
    # only drops rows when `temp` is in the formula.
    kept <- as.integer(rownames(fit$predict_data)[1:3])
    p_raw <- predict(fit, newdata = wine[kept, c("x", "temp", "contact")])
    p_fit <- predict(fit)
    expect_equal(unclass(p_raw)[, ], p_fit[1:3, ], tolerance = 1e-10, ignore_attr = TRUE)
  }
})
