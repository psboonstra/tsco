test_that("joint tests work for a single-predictor model", {
  dat <- make_wine_test_data()

  fit <- tsco(
    rating ~ temp,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  wald <- summary(fit, joint_test = "Wald")$joint_tests

  expect_identical(rownames(lrt), "temp")
  expect_true(is.finite(lrt[["Chisq"]]))
  expect_true(is.finite(lrt[["Pr(>Chisq)"]]))
  expect_true(lrt[["Chisq"]] >= 0)

  # PO | C | MR with K = 5 and C at the fourth level gives K - C + 1 = 3
  # degrees of freedom, matching the manuscript.
  expect_equal(lrt[["df"]], 3)
  expect_equal(wald[["df"]], 3L)
})

test_that("the intercept-only reduced fit matches a manual null log-likelihood", {
  dat <- make_wine_test_data()

  fit <- tsco(
    rating ~ temp,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  # An intercept-only PO model is saturated in the category probabilities, so
  # its maximum-likelihood fit is the observed category distribution.
  null_loglik <- function(counts) {
    p <- counts / sum(counts)
    sum(counts[counts > 0] * log(p[counts > 0]))
  }

  y <- as.character(dat$individual$rating)
  lower <- fit$lower_levels

  counts1 <- table(factor(y[y %in% lower], levels = lower))
  counts2 <- table(
    factor(
      ifelse(y %in% lower, fit$lower_collapsed_label, y),
      levels = fit$collapsed_levels
    )
  )

  expected <- 2 * (
    as.numeric(logLik(fit)) -
      (null_loglik(as.numeric(counts1)) + null_loglik(as.numeric(counts2)))
  )

  expect_equal(
    summary(fit, joint_test = "LRT")$joint_tests[["Chisq"]],
    expected,
    tolerance = 1e-6
  )
})

test_that("case weights work under the default engine", {
  dat <- make_wine_test_data()

  args <- list(
    formula = rating ~ temp + contact,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  fit <- do.call(tsco, args)
  fit_w <- do.call(
    tsco,
    c(args, list(weights = rep(2, nrow(dat$individual))))
  )

  # Doubling every weight doubles the log-likelihood and leaves the estimates
  # unchanged, whichever backend the weighted fit ends up using.
  expect_equal(
    as.numeric(logLik(fit_w)),
    2 * as.numeric(logLik(fit)),
    tolerance = 1e-4
  )

  expect_equal(
    unname(coef(fit_w)),
    unname(coef(fit)),
    tolerance = 1e-4
  )
})

test_that("reported coefficient names do not depend on the fitting engine", {
  dat <- make_wine_test_data()

  args <- list(
    formula = rating ~ temp + contact,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  fit_orm <- do.call(tsco, c(args, list(po_engine = "orm")))
  fit_vglm <- do.call(tsco, c(args, list(po_engine = "vglm")))

  expect_identical(fit_orm$stage1_engine, "orm")
  expect_identical(fit_vglm$stage1_engine, "vglm")

  expect_identical(names(coef(fit_orm)), names(coef(fit_vglm)))
  expect_identical(rownames(vcov(fit_orm)), rownames(vcov(fit_vglm)))
  expect_equal(unname(coef(fit_orm)), unname(coef(fit_vglm)), tolerance = 1e-4)

  # Thresholds are reported on the lower-tail scale, so they must not carry
  # orm's upper-tail `y>=` labels.
  expect_false(any(grepl(">=", names(coef(fit_orm)), fixed = TRUE)))

  expect_identical(
    names(coef(fit_orm, stage = "stage1"))[seq_along(fit_orm$lower_levels[-1])],
    paste0("Y<=", utils::head(fit_orm$lower_levels, -1L))
  )

  # The stage-2 collapsed category is reported as `Y < C`, not by its
  # internal placeholder label.
  expect_identical(
    names(coef(fit_orm, stage = "stage2"))[1L],
    paste0("Y<", fit_orm$cutoff_level)
  )
})

test_that("inline formula transformations are supported", {
  dat <- make_wine_test_data()

  wine <- dat$individual
  set.seed(408)
  wine$z <- stats::runif(nrow(wine), 1, 10)
  wine$log_z <- log(wine$z)

  fit_inline <- tsco(
    rating ~ log(z) + temp,
    data = wine,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  fit_precomputed <- tsco(
    rating ~ log_z + temp,
    data = wine,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_equal(
    as.numeric(logLik(fit_inline)),
    as.numeric(logLik(fit_precomputed))
  )

  expect_equal(
    unname(coef(fit_inline)),
    unname(coef(fit_precomputed)),
    tolerance = 1e-8
  )

  # `newdata` supplies the raw variable, not the transformed column.
  newdata <- data.frame(
    z = c(2, 8),
    temp = factor(c("cold", "warm"), levels = levels(wine$temp))
  )

  p_inline <- predict(fit_inline, newdata = newdata)

  expect_equal(unname(rowSums(p_inline)), rep(1, 2), tolerance = 1e-10)

  expect_equal(
    unname(p_inline),
    unname(predict(
      fit_precomputed,
      newdata = transform(newdata, log_z = log(newdata$z))
    )),
    tolerance = 1e-8
  )

  expect_true(
    all(is.finite(summary(fit_inline, joint_test = "LRT")$joint_tests$Chisq))
  )
})

test_that("a plain formula is unaffected by right-hand-side reconstruction", {
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

  expect_identical(
    attr(stats::terms(fit$stage1_formula), "term.labels"),
    c("temp", "contact")
  )

  expect_identical(
    attr(stats::terms(fit$stage2_formula), "term.labels"),
    c("temp", "contact")
  )
})

test_that("interactions survive right-hand-side reconstruction", {
  dat <- make_wine_test_data()

  fit <- tsco(
    rating ~ temp * contact,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_identical(
    attr(stats::terms(fit$stage1_formula), "term.labels"),
    c("temp", "contact", "temp:contact")
  )

  expect_equal(unname(rowSums(predict(fit))), rep(1, nrow(dat$individual)),
               tolerance = 1e-10)
})

test_that("separation is reported in the summary diagnostics", {
  dat <- make_wine_test_data()

  # Rating 5 never occurs under cold temperature, which quasi-separates the
  # multinomial stages.
  fit <- tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "multinomial",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  s <- summary(fit, joint_test = "none")

  expect_true(s$diagnostics$stage1$separated ||
                s$diagnostics$stage2$separated)

  expect_output(print(s), "Fit diagnostics")
  expect_output(print(s), "separation")

  # A well-behaved fit reports nothing.
  clean <- tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  s_clean <- summary(clean, joint_test = "none")

  expect_false(s_clean$diagnostics$stage1$separated)
  expect_false(s_clean$diagnostics$stage2$separated)
  expect_true(s_clean$diagnostics$stage1$converged)
})

test_that("printed coefficient tables state the reported scale", {
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

  out <- utils::capture.output(print(summary(fit, joint_test = "none")))

  expect_true(any(grepl("alpha_k - x'beta", out, fixed = TRUE)))
  expect_true(any(grepl("toward lower categories", out, fixed = TRUE)))
  expect_true(any(grepl("component BICs do not sum", out, fixed = TRUE)))
})

test_that("stage fits do not store a copy of the analysis data in their call", {
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

  call1 <- fit$fit_stage1$call
  call2 <- fit$fit_stage2@call

  expect_true(is.name(call1$data))
  expect_true(is.name(call2$data))
  expect_lt(as.numeric(utils::object.size(call1)), 1e5)
  expect_lt(as.numeric(utils::object.size(call2)), 1e5)
})

test_that("predicted probability columns are matched by name, not position", {
  dat <- make_wine_test_data()

  fit <- tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  # `.tsco_align_prob()` warns whenever it has to fall back to column order;
  # neither backend should trigger that for a well-formed fit.
  expect_silent(predict(fit))

  p <- predict(fit)
  expect_identical(colnames(p), fit$levels)
})
