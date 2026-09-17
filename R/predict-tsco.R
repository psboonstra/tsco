#' Predict from a TSCO model
#'
#' Computes predicted probabilities or predicted classes from a fitted
#' two-stage conditional odds model.
#'
#' For lower-partition categories, unconditional probabilities are computed as
#'
#' \deqn{
#'   Pr(Y = k | X) =
#'   Pr(Y < C | X) Pr(Y = k | Y < C, X).
#' }
#'
#' For upper-partition categories, probabilities are taken directly from the
#' stage-2 collapsed model.
#'
#' If a prediction row contains a factor level that was observed in the full
#' analysis data but not represented among the observations contributing to the
#' stage-1 model, then the conditional lower-partition probabilities are not
#' identifiable for that row. The `unsupported_stage1` argument controls how
#' such rows are handled.
#'
#' Stage-1 support is currently checked marginally for each factor: a row is
#' flagged when one of its factor levels was absent from the stage-1 fitting
#' data. This check does not detect an unobserved combination of otherwise
#' observed factor levels in a model containing interactions such as `x1 * x2`.
#' Predictions for such combinations may therefore be treated as supported and
#' should be regarded as model-based extrapolations.
#'
#' @param object An object of class `"tsco"`.
#' @param newdata Optional data frame containing predictor values at which to
#'   predict. If omitted, predictions are returned for the original observations
#'   or grouped covariate patterns.
#' @param type Prediction type. `"prob"` and `"response"` return unconditional
#'   category probabilities. `"class"` returns the most probable outcome
#'   category. `"stage"` returns stage-specific probabilities, combined
#'   unconditional probabilities, the stage-2 lower-partition mass, and an
#'   indicator of rows unsupported by the stage-1 model.
#' @param unsupported_stage1 How to handle prediction rows whose factor levels
#'   were not represented in the stage-1 fitting data. Options are:
#'   `"partial"`: return `NA` for unidentified lower-category probabilities
#'   but retain identifiable upper-category probabilities;
#'   `"uniform"`: split the stage-2 lower-partition probability equally across
#'   lower categories;
#'   `"empirical"`: split the lower-partition probability according to the
#'   empirical lower-category distribution among stage-1 observations;
#'   `"na"`: return `NA` for the entire probability row.
#' @param ... Must be empty. Backend-specific prediction arguments are not
#'   supported because a TSCO model can use different engines across stages.
#'
#' @return If `type = "prob"` or `"response"`, a matrix of predicted
#'   probabilities with one column per original outcome level. If
#'   `type = "class"`, an ordered factor. If `type = "stage"`, a list with
#'   stage-specific probabilities, combined probabilities, lower-partition
#'   mass, `unsupported_stage1`, and the unsupported-row handling method.
#'
#' @method predict tsco
#' @export
predict.tsco <- function(
    object,
    newdata = NULL,
    type = c("prob", "response", "class", "stage"),
    unsupported_stage1 = c("partial", "uniform", "empirical", "na"),
    ...) {

  type <- match.arg(type)
  unsupported_stage1 <- match.arg(unsupported_stage1)
  dots <- list(...)

  if (length(dots) > 0L) {
    stop(
      "Additional prediction arguments are not supported; `...` must be empty.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------
  # Prepare prediction data.
  # --------------------------------------------------------------------------

  pred_data <- .tsco_prepare_newdata(object, newdata)

  unsupported <- .tsco_stage1_unsupported_rows(object, pred_data)
  supported <- !unsupported

  if (any(unsupported)) {
    msg <- switch(
      unsupported_stage1,
      partial = paste0(
        "Some prediction rows contain factor levels not represented in the ",
        "stage-1 fitting data. Conditional lower-partition probabilities are ",
        "not identifiable for these rows; lower-category probabilities are ",
        "returned as NA, while identifiable upper-category probabilities are ",
        "retained."
      ),
      uniform = paste0(
        "Some prediction rows contain factor levels not represented in the ",
        "stage-1 fitting data. Conditional lower-partition probabilities are ",
        "not identifiable for these rows; the lower-partition probability is ",
        "being split uniformly across lower categories."
      ),
      empirical = paste0(
        "Some prediction rows contain factor levels not represented in the ",
        "stage-1 fitting data. Conditional lower-partition probabilities are ",
        "not identifiable for these rows; the lower-partition probability is ",
        "being split according to the empirical lower-category distribution ",
        "among stage-1 observations."
      ),
      na = paste0(
        "Some prediction rows contain factor levels not represented in the ",
        "stage-1 fitting data. Conditional lower-partition probabilities are ",
        "not identifiable for these rows; combined probabilities are returned ",
        "as NA for those rows."
      )
    )

    warning(msg, call. = FALSE)
  }

  predict_stage_prob <- function(fit, engine, pred_data, levels) {
    .tsco_predict_stage_prob(
      fit,
      engine = engine,
      newdata = pred_data,
      levels = levels
    )
  }

  # --------------------------------------------------------------------------
  # Stage 2 predictions: available for all rows.
  # --------------------------------------------------------------------------

  p2 <- predict_stage_prob(
    object$fit_stage2,
    engine = object$stage2_engine,
    pred_data = pred_data,
    levels = .tsco_stage_levels(object, 2L)
  )
  p2 <- .tsco_align_prob(p2, .tsco_stage_levels(object, 2L))
  p2 <- .tsco_expand_prob(
    p2, .tsco_stage_levels(object, 2L), object$collapsed_levels
  )

  n <- nrow(pred_data)

  if (nrow(p2) != n) {
    stop(
      "Stage 2 prediction returned an unexpected number of rows.",
      call. = FALSE
    )
  }

  lower_mass <- as.numeric(p2[, object$lower_collapsed_label])

  # --------------------------------------------------------------------------
  # Stage 1 predictions.
  # --------------------------------------------------------------------------

  p1 <- matrix(
    NA_real_,
    nrow = n,
    ncol = length(object$lower_levels),
    dimnames = list(rownames(pred_data), object$lower_levels)
  )

  # Predict stage 1 only for rows whose factor levels are represented in the
  # stage-1 fitting data.
  if (any(supported)) {
    p1_supported <- predict_stage_prob(
      object$fit_stage1,
      engine = object$stage1_engine,
      pred_data = pred_data[supported, , drop = FALSE],
      levels = .tsco_stage_levels(object, 1L)
    )

    p1_supported <- .tsco_align_prob(
      p1_supported,
      .tsco_stage_levels(object, 1L)
    )

    p1_supported <- .tsco_expand_prob(
      p1_supported, .tsco_stage_levels(object, 1L), object$lower_levels
    )

    if (nrow(p1_supported) != sum(supported)) {
      stop(
        "Stage 1 prediction returned an unexpected number of rows.",
        call. = FALSE
      )
    }

    p1[supported, ] <- p1_supported
  }

  # For unsupported rows, optionally complete the unidentified conditional
  # lower-partition probabilities.
  if (any(unsupported)) {
    if (unsupported_stage1 == "uniform") {
      p1[unsupported, ] <- matrix(
        1 / length(object$lower_levels),
        nrow = sum(unsupported),
        ncol = length(object$lower_levels),
        dimnames = list(
          rownames(pred_data)[unsupported],
          object$lower_levels
        )
      )
    }

    if (unsupported_stage1 == "empirical") {
      emp <- .tsco_empirical_lower_probs(object)

      p1[unsupported, ] <- matrix(
        rep(emp, each = sum(unsupported)),
        nrow = sum(unsupported),
        ncol = length(object$lower_levels),
        byrow = FALSE,
        dimnames = list(
          rownames(pred_data)[unsupported],
          object$lower_levels
        )
      )
    }
  }

  # --------------------------------------------------------------------------
  # Combine stage probabilities.
  # --------------------------------------------------------------------------

  K <- length(object$levels)

  prob <- matrix(
    NA_real_,
    nrow = n,
    ncol = K,
    dimnames = list(rownames(pred_data), object$levels)
  )

  completed <- supported |
    (unsupported & unsupported_stage1 %in% c("uniform", "empirical"))

  partial <- unsupported & unsupported_stage1 == "partial"

  # Rows with identifiable or completed lower-partition conditional
  # probabilities get full probability vectors.
  if (any(completed)) {
    prob[completed, object$lower_levels] <-
      p1[completed, object$lower_levels, drop = FALSE] *
      as.numeric(lower_mass[completed])

    prob[completed, object$upper_levels] <-
      p2[completed, object$upper_levels, drop = FALSE]
  }

  # Rows with partial predictions retain identifiable upper probabilities.
  # Lower categories remain NA because only their sum is identified.
  if (any(partial)) {
    prob[partial, object$upper_levels] <-
      p2[partial, object$upper_levels, drop = FALSE]
  }

  # Rows with unsupported_stage1 == "na" remain all NA.

  # Numerical cleanup and normalization only for rows with complete probability
  # vectors.
  if (any(completed)) {
    tol <- sqrt(.Machine$double.eps)

    finite_idx <- is.finite(prob)

    prob[finite_idx & prob < 0 & prob > -tol] <- 0
    prob[finite_idx & prob > 1 & prob < 1 + tol] <- 1

    rs <- rowSums(prob[completed, , drop = FALSE])

    if (any(!is.finite(rs)) || any(rs <= 0)) {
      stop(
        "Predicted probabilities could not be normalized.",
        call. = FALSE
      )
    }

    prob[completed, ] <-
      prob[completed, , drop = FALSE] / rs
  }

  attr(prob, "unsupported_stage1") <- unsupported
  attr(prob, "unsupported_stage1_method") <- unsupported_stage1
  attr(prob, "lower_mass") <- lower_mass

  # --------------------------------------------------------------------------
  # Return requested output type.
  # --------------------------------------------------------------------------

  if (type %in% c("prob", "response")) {
    return(prob)
  }

  if (type == "class") {
    cls <- rep(NA_character_, n)

    if (unsupported_stage1 %in% c("uniform", "empirical", "na")) {
      ok <- stats::complete.cases(prob)

      if (any(ok)) {
        cls[ok] <- object$levels[
          max.col(prob[ok, , drop = FALSE], ties.method = "first")
        ]
      }
    }

    if (unsupported_stage1 == "partial") {
      # Supported rows have complete probability vectors.
      if (any(supported)) {
        cls[supported] <- object$levels[
          max.col(prob[supported, , drop = FALSE], ties.method = "first")
        ]
      }

      # For unsupported rows, class can still be identified if the largest
      # upper-category probability is larger than the entire unidentified
      # lower-partition mass. Otherwise, a lower category could be the class,
      # so leave as NA.
      if (any(unsupported)) {
        upper_prob <- p2[unsupported, object$upper_levels, drop = FALSE]
        max_upper <- apply(upper_prob, 1L, max)
        which_upper <- max.col(upper_prob, ties.method = "first")

        identifiable_upper <- max_upper > lower_mass[unsupported]

        unsupported_indices <- which(unsupported)

        cls[
          unsupported_indices[identifiable_upper]
        ] <- object$upper_levels[
          which_upper[identifiable_upper]
        ]
      }
    }

    return(
      factor(
        cls,
        levels = object$levels,
        ordered = TRUE
      )
    )
  }

  if (type == "stage") {
    return(
      list(
        stage1_conditional = p1,
        stage2_collapsed = p2,
        prob = prob,
        lower_mass = lower_mass,
        unsupported_stage1 = unsupported,
        unsupported_stage1_method = unsupported_stage1
      )
    )
  }
}
