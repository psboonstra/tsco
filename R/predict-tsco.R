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
#' @param object An object of class `"tsco"`.
#' @param newdata Optional data frame containing predictor values at which to
#'   predict. If omitted, predictions are returned for the original observations
#'   or grouped covariate patterns.
#' @param type Prediction type. `"prob"` and `"response"` return unconditional
#'   category probabilities. `"class"` returns the most probable outcome
#'   category. `"stage"` returns stage-specific probabilities and combined
#'   unconditional probabilities.
#' @param ... Additional arguments passed to `predict()` for the underlying
#'   `VGAM::vglm` fits.
#'
#' @return If `type = "prob"` or `"response"`, a matrix of predicted
#'   probabilities with one column per original outcome level. If
#'   `type = "class"`, an ordered factor. If `type = "stage"`, a list with
#'   stage-specific and combined probabilities.
#'
#' @method predict tsco
#' @export
predict.tsco <- function(
    object,
    newdata = NULL,
    type = c("prob", "response", "class", "stage"),
    ...) {

  type <- match.arg(type)

  # If newdata is not supplied, predict for the original full analysis data.
  # Important: stage 1 was fit only among Y < C, so calling predict() on the
  # stage-1 fit without newdata would return predictions only for that subset.
  if (is.null(newdata)) {
    if (!is.null(object$predict_data)) {
      pred_data <- object$predict_data
    } else {
      pred_data <- NULL
    }
  } else {
    pred_data <- newdata
  }

  predict_stage_prob <- function(fit, pred_data, ...) {
    if (is.null(pred_data)) {
      VGAM::predictvglm(fit, type = "response", ...)
    } else {
      VGAM::predictvglm(fit, newdata = pred_data, type = "response", ...)
    }
  }

  p1 <- predict_stage_prob(object$fit_stage1, pred_data, ...)
  p2 <- predict_stage_prob(object$fit_stage2, pred_data, ...)

  p1 <- .tsco_align_prob(p1, object$lower_levels)
  p2 <- .tsco_align_prob(p2, object$collapsed_levels)

  if (nrow(p1) != nrow(p2)) {
    stop(
      "Stage 1 and stage 2 predictions have different numbers of rows. ",
      "If this object was created with an older version of tsco(), refit the ",
      "model or supply `newdata` explicitly.",
      call. = FALSE
    )
  }

  n <- nrow(p2)
  K <- length(object$levels)

  prob <- matrix(
    0,
    nrow = n,
    ncol = K,
    dimnames = list(rownames(p2), object$levels)
  )

  p_lower <- p2[, object$lower_collapsed_label]

  # Lower categories:
  # Pr(Y = k | X) = Pr(Y < C | X) Pr(Y = k | Y < C, X)
  prob[, object$lower_levels] <-
    p1[, object$lower_levels, drop = FALSE] * as.numeric(p_lower)

  # Upper categories:
  # These are modeled directly in the collapsed stage-2 model.
  prob[, object$upper_levels] <-
    p2[, object$upper_levels, drop = FALSE]

  # Numerical cleanup.
  tol <- sqrt(.Machine$double.eps)

  prob[prob < 0 & prob > -tol] <- 0
  prob[prob > 1 & prob < 1 + tol] <- 1

  rs <- rowSums(prob)

  if (any(!is.finite(rs)) || any(rs <= 0)) {
    stop(
      "Predicted probabilities could not be normalized.",
      call. = FALSE
    )
  }

  prob <- prob / rs

  if (type %in% c("prob", "response")) {
    return(prob)
  }

  if (type == "class") {
    cls <- object$levels[max.col(prob, ties.method = "first")]

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
        prob = prob
      )
    )
  }
}
