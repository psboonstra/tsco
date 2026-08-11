# Benchmark proportional-odds backends for a high-cardinality conditional stage.
#
# The simulated outcome represents discretized viral load below a clinical
# cutoff: many ordered lower categories and continuous covariates. The default
# 50- and 100-level settings are sized to run within typical workstation memory.
# Run interactively with:
# source("dev/benchmark-orm-vgam-high-cardinality.R")
#
# A 200-level stress run can be requested explicitly, but it may require much
# more memory for VGAM; increase n_per_level enough to avoid empty categories.

simulate_conditional_stage <- function(n_levels, n_per_level) {
  n <- n_levels * n_per_level
  x_age <- stats::rnorm(n)
  x_treatment <- stats::rbinom(n, size = 1L, prob = 0.5)
  linear_predictor <- 0.35 * x_age - 0.55 * x_treatment
  cutpoints <- seq(4, -4, length.out = n_levels - 1L)

  p_ge <- vapply(
    cutpoints,
    function(cutpoint) stats::plogis(cutpoint - linear_predictor),
    numeric(n)
  )

  probabilities <- cbind(
    1 - p_ge[, 1],
    p_ge[, -ncol(p_ge), drop = FALSE] - p_ge[, -1, drop = FALSE],
    p_ge[, ncol(p_ge)]
  )

  category <- apply(
    probabilities,
    1,
    function(probability) sample.int(n_levels, size = 1L, prob = probability)
  )

  data.frame(
    x_age = x_age,
    x_treatment = x_treatment,
    y = ordered(category, levels = seq_len(n_levels))
  )
}

measure_fit <- function(fitter, data) {
  gc(reset = TRUE)
  start <- proc.time()[["elapsed"]]
  fit <- fitter(data)
  elapsed_seconds <- proc.time()[["elapsed"]] - start
  memory <- gc()

  list(
    fit = fit,
    elapsed_seconds = elapsed_seconds,
    peak_mb = memory["Vcells", "max used"] * 8 / 1024^2
  )
}

benchmark_conditional_stage <- function(
    n_levels = c(50L, 100L),
    n_per_level = 25L,
    repetitions = 1L,
    seed = 6187L) {
  stopifnot(
    all(n_levels >= 3L),
    n_per_level >= 1L,
    repetitions >= 1L
  )

  set.seed(seed)
  fitters <- list(
    orm = function(data) {
      rms::orm(y ~ x_age + x_treatment, data = data, family = "logistic")
    },
    vgam = function(data) {
      VGAM::vglm(
        y ~ x_age + x_treatment,
        family = VGAM::cumulative(
          link = "logitlink",
          parallel = TRUE,
          reverse = FALSE
        ),
        data = data
      )
    }
  )

  measurements <- list()
  agreement <- list()
  row <- 1L

  for (levels in n_levels) {
    data <- simulate_conditional_stage(levels, n_per_level)
    fits <- list()

    for (replicate in seq_len(repetitions)) {
      engines <- if (replicate %% 2L == 1L) c("orm", "vgam") else c("vgam", "orm")

      for (engine in engines) {
        result <- measure_fit(fitters[[engine]], data)
        fits[[engine]] <- result$fit
        measurements[[row]] <- data.frame(
          n_levels = levels,
          n_observations = nrow(data),
          repetition = replicate,
          engine = engine,
          elapsed_seconds = result$elapsed_seconds,
          peak_mb = result$peak_mb
        )
        row <- row + 1L
      }
    }

    orm_probability <- stats::predict(fits$orm, type = "fitted.ind")
    vgam_probability <- VGAM::fittedvlm(fits$vgam)
    agreement[[length(agreement) + 1L]] <- data.frame(
      n_levels = levels,
      log_likelihood_difference = abs(
        as.numeric(stats::logLik(fits$orm)) - fits$vgam@criterion$loglikelihood
      ),
      max_probability_difference = max(abs(orm_probability - vgam_probability))
    )
  }

  measurements <- do.call(rbind, measurements)
  summary <- do.call(
    rbind,
    lapply(split(measurements, list(measurements$n_levels, measurements$engine), drop = TRUE), function(x) {
      data.frame(
        n_levels = x$n_levels[[1]],
        n_observations = x$n_observations[[1]],
        engine = x$engine[[1]],
        median_elapsed_seconds = stats::median(x$elapsed_seconds),
        min_elapsed_seconds = min(x$elapsed_seconds),
        max_elapsed_seconds = max(x$elapsed_seconds),
        median_peak_mb = stats::median(x$peak_mb)
      )
    })
  )
  rownames(summary) <- NULL
  summary <- summary[order(summary$n_levels, summary$engine), ]

  list(
    measurements = measurements,
    summary = summary,
    agreement = do.call(rbind, agreement),
    settings = list(
      n_levels = n_levels,
      n_per_level = n_per_level,
      repetitions = repetitions,
      seed = seed,
      R = R.version.string,
      rms = as.character(utils::packageVersion("rms")),
      VGAM = as.character(utils::packageVersion("VGAM"))
    )
  )
}

benchmark_results <- benchmark_conditional_stage()
benchmark_results$summary
benchmark_results$agreement
