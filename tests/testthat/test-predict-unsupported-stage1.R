make_sparse_stage1_data <- function() {
  grouped <- data.frame(
    grp = factor(c("A", "B", "U"), levels = c("A", "B", "U")),
    a = c(18, 8, 0),
    b = c(7, 17, 0),
    c = c(8, 6, 18),
    d = c(4, 9, 12)
  )

  counts <- as.matrix(grouped[c("a", "b", "c", "d")])
  response_levels <- colnames(counts)

  individual <- data.frame(
    grp = factor(
      rep(grouped$grp, times = rowSums(counts)),
      levels = levels(grouped$grp)
    ),
    y = ordered(
      unlist(
        lapply(
          seq_len(nrow(counts)),
          function(i) rep(response_levels, times = counts[i, ])
        ),
        use.names = FALSE
      ),
      levels = response_levels
    )
  )

  list(grouped = grouped, individual = individual)
}


fit_sparse_stage1_models <- function() {
  dat <- make_sparse_stage1_data()

  list(
    grouped = tsco(
      cbind(a, b, c, d) ~ grp,
      data = dat$grouped,
      cutoff_level = "c",
      stage1 = "po",
      stage2 = "po",
      po_engine = "vglm",
      warn_degenerate = FALSE
    ),
    individual = tsco(
      y ~ grp,
      data = dat$individual,
      cutoff_level = "c",
      stage1 = "po",
      stage2 = "po",
      po_engine = "vglm",
      warn_degenerate = FALSE
    )
  )
}


predict_sparse_modes <- function(fit, newdata) {
  methods <- c("partial", "uniform", "empirical", "na")

  stats::setNames(
    lapply(
      methods,
      function(method) {
        suppressWarnings(
          predict(
            fit,
            newdata = newdata,
            type = "stage",
            unsupported_stage1 = method
          )
        )
      }
    ),
    methods
  )
}


test_that("default orm handles a stage-1 factor level with no observations", {
  skip_if_not_installed("rms")

  dat <- make_sparse_stage1_data()
  fit <- tsco(
    y ~ grp,
    data = dat$individual,
    cutoff_level = "c",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_identical(fit$stage1_engine, "orm")
  expect_s3_class(fit$fit_stage1, "orm")
  expect_identical(levels(fit$data_stage1$grp), c("A", "B"))
  expect_identical(fit$stage1_factor_levels$grp, c("A", "B"))

  newdata <- data.frame(
    grp = factor(c("A", "U"), levels = c("A", "B", "U"))
  )
  modes <- predict_sparse_modes(fit, newdata)

  for (method in names(modes)) {
    expect_equal(dim(modes[[method]]$prob), c(2L, 4L))
    expect_identical(modes[[method]]$unsupported_stage1, c(FALSE, TRUE))
  }

  expect_equal(unname(rowSums(modes$uniform$prob)), c(1, 1), tolerance = 1e-10)
  expect_equal(unname(rowSums(modes$empirical$prob)), c(1, 1), tolerance = 1e-10)
  expect_true(all(is.na(modes$partial$prob[2L, fit$lower_levels])))
  expect_true(all(is.na(modes$na$prob[2L, ])))

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  wald <- summary(fit, joint_test = "Wald")$joint_tests
  expect_identical(rownames(lrt), "grp")
  expect_identical(rownames(wald), "grp")
  expect_equal(sum(!is.finite(lrt$Chisq)), 0L)
  expect_equal(sum(!is.finite(wald$Chisq)), 0L)
})


test_that("unsupported stage-1 modes work for grouped and individual data", {
  skip_if_not_installed("VGAM")

  fits <- fit_sparse_stage1_models()
  newdata <- data.frame(
    grp = factor(c("A", "U"), levels = c("A", "B", "U"))
  )
  results <- lapply(fits, predict_sparse_modes, newdata = newdata)
  expected_empirical <- c(a = 26, b = 24) / 50

  for (representation in names(fits)) {
    fit <- fits[[representation]]
    modes <- results[[representation]]
    supported <- 1L
    unsupported <- 2L

    expect_identical(modes$partial$unsupported_stage1, c(FALSE, TRUE))
    expect_identical(
      attr(modes$partial$prob, "unsupported_stage1"),
      c(FALSE, TRUE)
    )
    expect_identical(
      attr(modes$partial$prob, "unsupported_stage1_method"),
      "partial"
    )
    expect_equal(
      attr(modes$partial$prob, "lower_mass"),
      modes$partial$lower_mass
    )
    expect_true(all(is.na(
      modes$partial$prob[unsupported, fit$lower_levels]
    )))
    expect_equal(
      modes$partial$prob[unsupported, fit$upper_levels],
      modes$partial$stage2_collapsed[unsupported, fit$upper_levels]
    )

    expect_equal(
      unname(modes$uniform$prob[unsupported, fit$lower_levels]),
      rep(
        modes$uniform$lower_mass[unsupported] / length(fit$lower_levels),
        length(fit$lower_levels)
      )
    )
    expect_equal(
      unname(modes$empirical$prob[unsupported, fit$lower_levels]),
      unname(modes$empirical$lower_mass[unsupported] * expected_empirical)
    )
    expect_equal(
      unname(rowSums(modes$uniform$prob)),
      rep(1, nrow(newdata)),
      tolerance = 1e-10
    )
    expect_equal(
      unname(rowSums(modes$empirical$prob)),
      rep(1, nrow(newdata)),
      tolerance = 1e-10
    )

    expect_true(all(is.na(modes$na$prob[unsupported, ])))

    for (method in names(modes)) {
      expect_equal(
        modes[[method]]$prob[supported, ],
        modes$partial$prob[supported, ]
      )
    }
  }

  for (method in names(results$grouped)) {
    expect_equal(
      results$grouped[[method]]$prob,
      results$individual[[method]]$prob,
      tolerance = 1e-5,
      ignore_attr = TRUE
    )
  }
})


test_that("partial classes require an upper probability above lower mass", {
  skip_if_not_installed("VGAM")

  fits <- fit_sparse_stage1_models()
  newdata <- data.frame(
    grp = factor(c("B", "U"), levels = c("A", "B", "U"))
  )

  for (representation in names(fits)) {
    fit <- fits[[representation]]

    # Treat both rows as unsupported to isolate both branches of the class rule.
    fit$stage1_factor_levels$grp <- "A"

    stage <- suppressWarnings(
      predict(
        fit,
        newdata = newdata,
        type = "stage",
        unsupported_stage1 = "partial"
      )
    )
    class <- suppressWarnings(
      predict(
        fit,
        newdata = newdata,
        type = "class",
        unsupported_stage1 = "partial"
      )
    )

    upper_prob <- stage$stage2_collapsed[, fit$upper_levels, drop = FALSE]
    identifiable <- apply(upper_prob, 1L, max) > stage$lower_mass
    expected <- rep(NA_character_, nrow(newdata))
    expected[identifiable] <- fit$upper_levels[
      max.col(upper_prob[identifiable, , drop = FALSE], ties.method = "first")
    ]

    expect_identical(unname(identifiable), c(FALSE, TRUE))
    expect_identical(as.character(class), expected)
  }
})
