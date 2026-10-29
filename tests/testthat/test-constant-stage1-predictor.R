# A predictor that is constant among the observations below the cutoff -- the
# manuscript's own two-arm design when the treatment arm never falls below it
# -- is dropped from the stage-1 model rather than handed to the fitter as a
# one-level factor. Stage 2 keeps the full formula, the joint tests cover the
# term through stage 2 alone, and prediction treats the other values as
# unsupported in stage 1.

two_arm_grouped <- function() {
  data.frame(
    arm = factor(c("control", "treatment")),
    a = c(18, 0), b = c(14, 0), c = c(10, 18), d = c(8, 12)
  )
}

two_arm_individual <- function() {
  g <- two_arm_grouped()
  lev <- c("a", "b", "c", "d")
  data.frame(
    arm = factor(rep(g$arm, times = rowSums(g[, lev]))),
    y = ordered(
      unlist(lapply(seq_len(nrow(g)), function(i) rep(lev, times = unlist(g[i, lev])))),
      levels = lev
    )
  )
}

# By hand: stage 1 is the control arm's a:b split; stage 2 is a PO model on
# {<c, c, d} fitted with VGAM. See the test body for the combination.
hand_tsco <- function(g) {
  s1 <- c(a = g$a[1], b = g$b[1]) / (g$a[1] + g$b[1])
  s2 <- VGAM::vglm(
    cbind(a + b, c, d) ~ arm, data = g,
    family = VGAM::cumulative(link = "logitlink", parallel = TRUE)
  )
  p2 <- VGAM::predictvglm(s2, newdata = data.frame(arm = g$arm), type = "response")
  out <- cbind(a = p2[, 1] * s1["a"], b = p2[, 1] * s1["b"], c = p2[, 2], d = p2[, 3])
  dimnames(out) <- list(as.character(g$arm), c("a", "b", "c", "d"))
  list(prob = out, lower_mass = unname(p2[, 1]), s1 = s1, stage2 = s2)
}

fit_two_arm <- function(data, formula, cutoff_level = "c", ...) {
  tsco(formula, data = data, cutoff_level = cutoff_level, stage1 = "po",
       stage2 = "po", warn_degenerate = FALSE, ...)
}


test_that("a predictor constant below the cutoff is dropped from stage 1 with a warning", {
  g <- two_arm_grouped()

  expect_warning(
    fit <- fit_two_arm(g, cbind(a, b, c, d) ~ arm),
    "`arm` take a single value among the observations below the cutoff"
  )

  expect_identical(fit$stage1_dropped_terms, "arm")
  expect_identical(fit$term_labels, "arm")
  expect_length(attr(stats::terms(fit$stage1_formula), "term.labels"), 0L)
  expect_identical(attr(stats::terms(fit$stage2_formula), "term.labels"), "arm")
  expect_identical(fit$stage1_factor_levels$arm, "control")

  expect_identical(
    names(coef(fit)),
    c("stage1:Y<=a", "stage2:Y<c", "stage2:Y<=c", "stage2:armtreatment")
  )

  # Stage 1 is intercept-only and saturated: its threshold is logit(18/32).
  expect_equal(unname(coef(fit)["stage1:Y<=a"]), qlogis(18 / 32), tolerance = 1e-6)
})


test_that("an intercept-only PO stage is fitted with VGAM whatever po_engine says", {
  ind <- two_arm_individual()

  fit_auto <- suppressWarnings(fit_two_arm(ind, y ~ arm))
  expect_identical(fit_auto$stage1_engine, "vglm")
  expect_identical(fit_auto$stage2_engine, "orm")

  fit_orm <- suppressWarnings(fit_two_arm(ind, y ~ arm, po_engine = "orm"))
  expect_identical(fit_orm$stage1_engine, "vglm")
  expect_identical(fit_orm$stage2_engine, "orm")

  expect_equal(coef(fit_orm), coef(fit_auto), tolerance = 1e-6)
})


test_that("two-arm predictions match the hand computation under every policy", {
  g <- two_arm_grouped()
  hand <- hand_tsco(g)
  nd <- data.frame(arm = factor(c("control", "treatment"), levels = levels(g$arm)))

  for (data in list(g, two_arm_individual())) {
    formula <- if (is.ordered(data$y)) y ~ arm else cbind(a, b, c, d) ~ arm
    fit <- suppressWarnings(fit_two_arm(data, formula))

    # Agreement is to solver tolerance: the hand fit is VGAM, the individual
    # stage 2 is orm.
    # Stage-2 mass below the cutoff is identified for both arms.
    expect_warning(p <- predict(fit, newdata = nd), "not represented in the stage-1")
    expect_equal(unname(attr(p, "lower_mass")), hand$lower_mass, tolerance = 1e-4)

    # partial: control fully predicted, treatment's lower split NA.
    expect_equal(unclass(p)["1", ], hand$prob["control", ], tolerance = 1e-4, ignore_attr = TRUE)
    expect_true(all(is.na(p["2", c("a", "b")])))
    expect_equal(unclass(p)["2", c("c", "d")], hand$prob["treatment", c("c", "d")],
                 tolerance = 1e-4, ignore_attr = TRUE)
    expect_identical(unname(attr(p, "unsupported_stage1")), c(FALSE, TRUE))

    # empirical: borrows the control arm's a:b split, which is the hand value.
    pe <- suppressWarnings(predict(fit, newdata = nd, unsupported_stage1 = "empirical"))
    expect_equal(unclass(pe)[, ], hand$prob, tolerance = 1e-4, ignore_attr = TRUE)

    # uniform: equal split of the same mass.
    pu <- suppressWarnings(predict(fit, newdata = nd, unsupported_stage1 = "uniform"))
    expect_equal(unname(pu["2", c("a", "b")]), rep(hand$lower_mass[2] / 2, 2), tolerance = 1e-4)

    # na: whole treatment row NA.
    pn <- suppressWarnings(predict(fit, newdata = nd, unsupported_stage1 = "na"))
    expect_true(all(is.na(pn["2", ])))
    expect_false(anyNA(pn["1", ]))

    # Class under "partial": treatment's top category beats the whole lower mass.
    cls <- suppressWarnings(predict(fit, newdata = nd, type = "class"))
    expect_identical(as.character(cls), c("a", "d"))
  }
})


test_that("the joint tests for a term missing from stage 1 come from stage 2 alone", {
  g <- two_arm_grouped()
  fit <- suppressWarnings(fit_two_arm(g, cbind(a, b, c, d) ~ arm))
  hand <- hand_tsco(g)

  counts <- cbind(g$a + g$b, g$c, g$d)
  ll <- function(f) sum(counts * log(VGAM::predictvglm(f, type = "response")))
  s2_null <- VGAM::vglm(
    cbind(a + b, c, d) ~ 1, data = g,
    family = VGAM::cumulative(link = "logitlink", parallel = TRUE)
  )
  manual_lrt <- 2 * (ll(hand$stage2) - ll(s2_null))

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  expect_identical(rownames(lrt), "arm")
  expect_equal(lrt$df, 1)
  expect_equal(lrt$Chisq, manual_lrt, tolerance = 1e-6)

  wald <- summary(fit, joint_test = "Wald")$joint_tests
  expect_identical(rownames(wald), "arm")
  expect_equal(wald$df, 1)
  b <- coef(fit)["stage2:armtreatment"]
  V <- vcov(fit)["stage2:armtreatment", "stage2:armtreatment"]
  expect_equal(wald$Chisq, unname(b^2 / V), tolerance = 1e-8)

  # The same statistics from the individual-level representation.
  fi <- suppressWarnings(fit_two_arm(two_arm_individual(), y ~ arm))
  expect_equal(summary(fi, joint_test = "LRT")$joint_tests$Chisq, manual_lrt, tolerance = 1e-6)
})


test_that("only the constant predictor's terms are dropped; others stay in stage 1", {
  ind <- two_arm_individual()
  set.seed(4)
  ind$w <- rnorm(nrow(ind))

  expect_warning(
    fit <- fit_two_arm(ind, y ~ arm * w),
    "term\\(s\\) `arm`, `arm:w` were dropped"
  )

  expect_identical(fit$stage1_dropped_terms, c("arm", "arm:w"))
  expect_identical(attr(stats::terms(fit$stage1_formula), "term.labels"), "w")
  expect_identical(fit$stage1_engine, "orm")
  expect_identical(
    names(coef(fit)),
    c("stage1:Y<=a", "stage1:w", "stage2:Y<c", "stage2:Y<=c",
      "stage2:armtreatment", "stage2:w", "stage2:armtreatment:w")
  )

  # The Wald table still covers every term of the user's formula, with
  # type-II term sets: arm = {arm, arm:w} in stage 2 only (2 df);
  # w = {w} in stage 1 plus {w, arm:w} in stage 2 (3 df); arm:w = 1 df.
  wald <- summary(fit, joint_test = "Wald")$joint_tests
  expect_identical(rownames(wald), c("arm", "w", "arm:w"))
  expect_equal(wald["arm", "df"], 2)
  expect_equal(wald["w", "df"], 3)
  expect_equal(wald["arm:w", "df"], 1)

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  expect_equal(lrt$df, c(2, 3, 1))
  expect_true(all(is.finite(lrt$Chisq)))

  # Prediction: control rows are supported, treatment rows are not.
  nd <- data.frame(arm = factor(c("control", "treatment")), w = c(0.5, 0.5))
  expect_warning(p <- predict(fit, newdata = nd), "not represented in the stage-1")
  expect_false(anyNA(p[1, ]))
  expect_true(all(is.na(p[2, c("a", "b")])))
})


test_that("a numeric predictor constant below the cutoff is dropped the same way", {
  ind <- two_arm_individual()
  ind$z <- ifelse(ind$arm == "control", 1, rnorm(nrow(ind)))

  expect_warning(fit <- fit_two_arm(ind, y ~ z), "`z` take a single value")
  expect_identical(fit$stage1_dropped_terms, "z")
  expect_length(attr(stats::terms(fit$stage1_formula), "term.labels"), 0L)
  expect_true("stage2:z" %in% names(coef(fit)))
})


test_that("multinomial stages handle a constant stage-1 predictor too", {
  ind <- two_arm_individual()
  fit <- suppressWarnings(tsco(y ~ arm, data = ind, cutoff_level = "c",
                               stage1 = "multinomial", stage2 = "multinomial",
                               warn_degenerate = FALSE))
  expect_identical(names(coef(fit))[1], "stage1:(Intercept):b")
  expect_equal(unname(coef(fit)[1]), log(14 / 18), tolerance = 1e-6)
  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  expect_identical(rownames(lrt), "arm")
  expect_equal(lrt$df, 2)
})


test_that("nothing changes when every predictor varies below the cutoff", {
  g <- two_arm_grouped()
  g$a[2] <- 5
  g$b[2] <- 4
  expect_no_warning(fit <- fit_two_arm(g, cbind(a, b, c, d) ~ arm))
  expect_length(fit$stage1_dropped_terms, 0L)
  expect_identical(attr(stats::terms(fit$stage1_formula), "term.labels"), "arm")
  expect_true("stage1:armtreatment" %in% names(coef(fit)))
})


# --- type-II joint tests ------------------------------------------------------

test_that("joint tests drop a term together with every term containing it", {
  ind <- two_arm_individual()
  set.seed(4)
  ind$w <- rnorm(nrow(ind))
  # Make stage 1 non-degenerate so every term is in both stages.
  ind$y[ind$arm == "treatment"][1:10] <- "a"

  fit <- function(f) tsco(f, data = ind, cutoff_level = "c", stage1 = "po",
                          stage2 = "po", warn_degenerate = FALSE)

  full <- fit(y ~ arm * w)
  lrt <- summary(full, joint_test = "LRT")$joint_tests

  # Manual nested fits on the shared constant-free scale.
  lr <- function(reduced) 2 * (full$logLik - reduced$logLik)
  ddf <- function(reduced) attr(logLik(full), "df") - attr(logLik(reduced), "df")

  no_arm <- fit(y ~ w)          # drops arm and arm:w
  no_w <- fit(y ~ arm)          # drops w and arm:w
  no_int <- fit(y ~ arm + w)    # drops arm:w only

  expect_equal(lrt["arm", "Chisq"], lr(no_arm), tolerance = 1e-6)
  expect_equal(lrt["arm", "df"], ddf(no_arm))
  expect_equal(lrt["w", "Chisq"], lr(no_w), tolerance = 1e-6)
  expect_equal(lrt["w", "df"], ddf(no_w))
  expect_equal(lrt["arm:w", "Chisq"], lr(no_int), tolerance = 1e-6)
  expect_equal(lrt["arm:w", "df"], ddf(no_int))
  expect_equal(lrt$df, c(4, 4, 2))

  wald <- summary(full, joint_test = "Wald")$joint_tests
  expect_equal(wald$df, c(4, 4, 2))

  # The non-hierarchical "drop the main effect only" model is never fitted,
  # so orm's marginality restriction no longer breaks the LRT.
  expect_true(all(is.finite(lrt$Chisq)))
})


test_that(".tsco_terms_containing follows variable-set containment", {
  f <- y ~ a * b + c + a:b:d
  expect_setequal(.tsco_terms_containing("a", f), c("a", "a:b", "a:b:d"))
  expect_setequal(.tsco_terms_containing("a:b", f), c("a:b", "a:b:d"))
  expect_identical(.tsco_terms_containing("c", f), "c")
  expect_identical(.tsco_terms_containing("d", f), character())   # not a term
  expect_identical(.tsco_terms_containing("a", y ~ c), character())
  expect_identical(.tsco_terms_containing("a", y ~ 1), character())
})


# --- the mirror image: no outcomes at or above the cutoff ---------------------

test_that("an arm with nothing at or above the cutoff is separation in stage 2, and is flagged", {
  g <- data.frame(arm = factor(c("control", "treatment")),
                  a = c(18, 20), b = c(14, 10), c = c(10, 0), d = c(8, 0))
  lev <- c("a", "b", "c", "d")
  ind <- data.frame(
    arm = factor(rep(g$arm, times = rowSums(g[, lev]))),
    y = ordered(unlist(lapply(1:2, function(i) rep(lev, times = unlist(g[i, lev])))), levels = lev)
  )

  for (engine in c("orm", "vglm")) {
    fit <- suppressWarnings(tsco(y ~ arm, data = ind, cutoff_level = "c",
                                 stage1 = "po", stage2 = "po", po_engine = engine,
                                 warn_degenerate = FALSE))

    # Nothing is dropped: the arm is present in stage 2, with an infinite effect.
    expect_length(fit$stage1_dropped_terms, 0L)
    expect_true("stage2:armtreatment" %in% names(coef(fit)))
    expect_gt(abs(coef(fit)["stage2:armtreatment"]), 10)

    s <- suppressWarnings(summary(fit, joint_test = "none"))
    expect_true(s$diagnostics$stage2$separated)
    expect_match(paste(s$diagnostics$stage2$messages, collapse = " "), "armtreatment")
    expect_false(s$diagnostics$stage1$separated)

    p <- predict(fit, newdata = data.frame(arm = factor(c("control", "treatment"))))
    # orm stops at a finite slope (~10), so the upper mass is ~1e-5, not 0.
    expect_lt(sum(p[2, c("c", "d")]), 1e-3)
    expect_equal(unname(p[2, c("a", "b")]), c(20, 10) / 30, tolerance = 1e-4)
    expect_false(anyNA(p))
  }
})


# --- nothing at all on one side of the cutoff is a clear error ----------------

test_that("a cutoff with no observations on one side is a clear error, with no stray warning", {
  g_hi <- data.frame(arm = factor(c("control", "treatment")),
                     a = c(18, 14), b = c(14, 12), c = c(10, 18), d = c(0, 0))
  g_lo <- data.frame(arm = factor(c("control", "treatment")),
                     a = c(0, 0), b = c(0, 0), c = c(10, 18), d = c(8, 6))
  lev <- c("a", "b", "c", "d")
  to_ind <- function(g) data.frame(
    arm = factor(rep(g$arm, times = rowSums(g[, lev]))),
    y = ordered(unlist(lapply(1:2, function(i) rep(lev, times = unlist(g[i, lev])))), levels = lev)
  )

  # No observation at or above the cutoff: stage 2 has one observed level.
  expect_no_warning(expect_error(
    fit_two_arm(g_hi, cbind(a, b, c, d) ~ arm, cutoff_level = "d"),
    "Stage 2 .* At least two are needed"
  ))
  expect_no_warning(expect_error(
    tsco(y ~ arm, data = to_ind(g_hi), levels = lev, cutoff_level = "d",
         stage1 = "po", stage2 = "po", warn_degenerate = FALSE),
    "Stage 2 .* At least two are needed"
  ))

  # No observation below the cutoff: stage 1 has no rows.
  expect_no_warning(expect_error(
    fit_two_arm(g_lo, cbind(a, b, c, d) ~ arm),
    "No grouped/count rows have positive counts below"
  ))
  expect_no_warning(expect_error(
    tsco(y ~ arm, data = to_ind(g_lo), levels = lev, cutoff_level = "c",
         stage1 = "po", stage2 = "po", warn_degenerate = FALSE),
    "No observations fall below `cutoff_level`"
  ))
})
