#' Fit a Two-Stage Conditional Odds (TSCO) Model
#'
#' Fits a two-stage conditional odds model for an ordinal outcome using
#' VGAM::vglm().
#'
#' The outcome is partitioned at a user-specified cutoff C.
#' Stage 1 models Y | Y < C. Stage 2 models the collapsed outcome
#' {Y < C, C, ..., K - 1}.
#'
#' The response may be either an individual-level ordinal outcome or a grouped
#' count response supplied using cbind(), as in
#' cbind(normal, mild, severe) ~ x.
#'
#' @param formula Model formula, e.g. y ~ x1 + x2 for individual-level data or
#'   cbind(y0, y1, y2) ~ x1 + x2 for grouped count data.
#' @param data A data.frame.
#' @param cutoff_level The first outcome level in the upper partition, supplied
#'   as a level label.
#' @param levels Optional character vector specifying the ordinal outcome level
#'   order. For individual-level data, observed outcome values must be contained
#'   in `levels`. For grouped count data with named cbind() columns, `levels`
#'   must contain exactly the response column names and is used to reorder the
#'   columns into ordinal order. For grouped count data without column names,
#'   `levels` is assumed to describe the existing cbind() column order.
#' @param stage1 Model for the conditional lower partition Y | Y < C.
#'   Either "po" for proportional odds or "multinomial".
#' @param stage2 Model for the collapsed marginal outcome \{Y < C, C, ..., K - 1\}.
#'   Either "multinomial" or "po".
#' @param po.reverse Logical. Passed to `VGAM::cumulative(reverse = ...)`.
#'   The default `FALSE` makes VGAM model lower cumulative probabilities,
#'   e.g. `Pr(Y < k | X)`, so that the fitted slope coefficients have the same
#'   sign as the manuscript parameterization
#'   `Pr(Y >= k | X) = expit(alpha_k - x' beta)`. If `TRUE`, VGAM models
#'   reverse cumulative probabilities and the slope coefficients are the
#'   negatives of the manuscript's beta parameters.
#' @param weights Optional numeric case weights.
#' @param na.action Missing-data action passed to model.frame().
#' @param warn_degenerate If TRUE, warn when K = 3, which reduces to two binary regressions.
#' @param ... Additional arguments passed to VGAM::vglm().
#'
#' @return An object of class "tsco".
#' @export
tsco <- function(
    formula,
    data,
    cutoff_level,
    levels = NULL,
    stage1 = c("po", "multinomial"),
    stage2 = c("multinomial", "po"),
    po.reverse = FALSE,
    weights = NULL,
    na.action = stats::na.omit,
    warn_degenerate = TRUE,
    ...) {

  if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop(
      "Package 'VGAM' is required. Install it with install.packages('VGAM').",
      call. = FALSE
    )
  }

  stage1 <- match.arg(stage1)
  stage2 <- match.arg(stage2)

  if (missing(data)) {
    data <- parent.frame()
  }

  if (length(formula) != 3L) {
    stop(
      "`formula` must be a two-sided formula, e.g. y ~ x1 + x2 or cbind(y0, y1, y2) ~ x.",
      call. = FALSE
    )
  }

  # Build model frame using original formula. This handles na.action and also
  # gives us a clean response object via model.response().
  mf <- stats::model.frame(
    formula = formula,
    data = data,
    na.action = na.action
  )

  y <- stats::model.response(mf)

  if (is.null(y)) {
    stop("Could not find response variable in `formula`.", call. = FALSE)
  }

  # Identify response column in model frame so we can construct stage-specific
  # data frames containing only predictors plus a new internal response.
  tt <- stats::terms(formula, data = data)
  response_col <- attr(tt, "response")

  if (is.null(response_col) || response_col < 1L) {
    stop("Could not identify the response in `formula`.", call. = FALSE)
  }

  cov_mf <- mf[, -response_col, drop = FALSE]

  # Detect grouped/count response.
  is_grouped <- is.matrix(y) || is.data.frame(y)

  if (is_grouped) {
    y_mat <- as.matrix(y)

    if (!is.numeric(y_mat)) {
      stop("Grouped/count response must be numeric.", call. = FALSE)
    }

    if (ncol(y_mat) < 2L) {
      stop("Grouped/count response must have at least two columns.", call. = FALSE)
    }

    if (any(!is.finite(y_mat))) {
      stop("Grouped/count response contains non-finite values.", call. = FALSE)
    }

    if (any(y_mat < 0)) {
      stop("Grouped/count response counts must be nonnegative.", call. = FALSE)
    }

    if (any(abs(y_mat - round(y_mat)) > sqrt(.Machine$double.eps))) {
      warning("Non-integer grouped/count response counts detected.", call. = FALSE)
    }

    if (any(rowSums(y_mat) <= 0)) {
      stop("Each grouped/count row must have positive total count.", call. = FALSE)
    }

    y_colnames <- colnames(y_mat)

    if (!is.null(levels)) {
      y_levels <- as.character(levels)

      if (length(y_levels) != ncol(y_mat)) {
        stop(
          "`levels` must have length equal to the number of response columns.",
          call. = FALSE
        )
      }

      if (!is.null(y_colnames) &&
          !anyNA(y_colnames) &&
          !any(y_colnames == "")) {

        y_colnames <- as.character(y_colnames)

        missing_levels <- setdiff(y_levels, y_colnames)
        extra_cols <- setdiff(y_colnames, y_levels)

        if (length(missing_levels) > 0L || length(extra_cols) > 0L) {
          stop(
            "`levels` must contain exactly the grouped response column names when ",
            "the cbind() response has column names.\n",
            "Missing from response: ",
            if (length(missing_levels) == 0L) "none" else paste(missing_levels, collapse = ", "),
            "\nExtra response columns: ",
            if (length(extra_cols) == 0L) "none" else paste(extra_cols, collapse = ", "),
            call. = FALSE
          )
        }

        # Important: reorder count matrix columns into the desired ordinal order.
        y_mat <- y_mat[, y_levels, drop = FALSE]

      } else {
        # No usable column names. In this case, levels are assumed to describe
        # the existing cbind() column order.
        colnames(y_mat) <- y_levels
      }

    } else {
      if (is.null(y_colnames) || anyNA(y_colnames) || any(y_colnames == "")) {
        stop(
          "For grouped/count responses, either the cbind() response must have ",
          "column names or `levels` must be supplied.",
          call. = FALSE
        )
      }

      y_levels <- as.character(y_colnames)
    }

    colnames(y_mat) <- y_levels

  } else {
    # Individual-level response.
    if (!is.null(levels)) {
      y_levels <- as.character(levels)
    } else if (is.ordered(y)) {
      y_levels <- base::levels(y)
    } else if (is.factor(y)) {
      warning(
        "Outcome is a factor but not ordered; using factor level order.",
        call. = FALSE
      )
      y_levels <- base::levels(y)
    } else if (is.numeric(y) || is.integer(y)) {
      y_levels <- as.character(sort(unique(y)))
    } else {
      stop(
        "For nonnumeric individual-level outcomes, please supply an ordered factor or explicit `levels`.",
        call. = FALSE
      )
    }

    y_chr <- as.character(y)

    unknown_y <- setdiff(unique(y_chr), y_levels)
    if (length(unknown_y) > 0L) {
      stop(
        "Observed outcome values not found in `levels`: ",
        paste(unknown_y, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (anyDuplicated(y_levels)) {
    stop("`levels` must not contain duplicate values.", call. = FALSE)
  }

  K <- length(y_levels)

  if (K < 3L) {
    stop(
      "The TSCO factorization requires at least 3 ordered outcome levels.",
      call. = FALSE
    )
  }

  if (K == 3L && isTRUE(warn_degenerate)) {
    warning(
      "With K = 3, TSCO reduces to two binary regressions; ",
      "the PO and multinomial choices are equivalent within each stage.",
      call. = FALSE
    )
  }

  # Locate cutoff.
  cutoff_index <- match(as.character(cutoff_level), y_levels)

  if (is.na(cutoff_index)) {
    stop(
      "`cutoff_level` must be one of the outcome levels: ",
      paste(y_levels, collapse = ", "),
      call. = FALSE
    )
  }

  lower_levels <- y_levels[seq_len(cutoff_index - 1L)]
  upper_levels <- y_levels[cutoff_index:K]

  if (length(lower_levels) < 2L) {
    stop(
      "Stage 1 requires at least two categories below `cutoff_level`.",
      call. = FALSE
    )
  }

  if (length(upper_levels) < 1L) {
    stop(
      "`cutoff_level` must leave at least one category in the upper partition.",
      call. = FALSE
    )
  }

  # Create a unique dummy label for the collapsed lower outcome.
  lower_collapsed_label <- paste0(".Y_lt_", make.names(y_levels[cutoff_index]))
  while (lower_collapsed_label %in% y_levels) {
    lower_collapsed_label <- paste0(".", lower_collapsed_label)
  }

  collapsed_levels <- c(lower_collapsed_label, upper_levels)

  # Validate and align weights.
  if (!is.null(weights)) {
    if (!is.numeric(weights)) {
      stop("`weights` must be a numeric vector.", call. = FALSE)
    }

    if (length(weights) == nrow(mf)) {
      w_all <- weights
    } else if (length(weights) == nrow(data)) {
      rn <- suppressWarnings(as.integer(rownames(mf)))

      if (anyNA(rn)) {
        stop(
          "Could not align `weights` with model frame. Ensure weight length equals nrow(model.frame).",
          call. = FALSE
        )
      }

      w_all <- weights[rn]
    } else {
      stop(
        "`weights` must equal nrow(data) or the number of complete cases.",
        call. = FALSE
      )
    }
  } else {
    w_all <- NULL
  }

  # Internal helper: construct formula with new response and original RHS.
  make_stage_formula <- function(response_name) {
    stats::as.formula(
      call("~", as.name(response_name), formula[[3L]]),
      env = environment(formula)
    )
  }

  # Helper to construct VGAM family.
  make_family <- function(kind) {
    if (kind == "po") {
      VGAM::cumulative(
        link = "logitlink",
        parallel = TRUE,
        reverse = po.reverse
      )
    } else if (kind == "multinomial") {
      VGAM::multinomial(refLevel = 1)
    }
  }

  # Model fitter.
  fit_vglm <- function(form, fam, dat, w, ...) {
    args <- list(
      formula = form,
      family = fam,
      data = dat
    )

    if (!is.null(w)) {
      args$weights <- w
    }

    args <- c(args, list(...))
    do.call(VGAM::vglm, args)
  }

  # --------------------------------------------------------------------------
  # Stage-specific data construction.
  # --------------------------------------------------------------------------

  if (is_grouped) {
    lower_idx <- match(lower_levels, y_levels)
    upper_idx <- match(upper_levels, y_levels)

    lower_counts <- y_mat[, lower_idx, drop = FALSE]
    upper_counts <- y_mat[, upper_idx, drop = FALSE]

    lower_total <- rowSums(lower_counts)

    n_groups <- nrow(y_mat)
    n_obs <- sum(rowSums(y_mat))
    n_stage1_groups <- sum(lower_total > 0)
    n_stage1_obs <- sum(lower_total)

    # Stage 1: conditional model among Y < C.
    # Rows with zero lower counts contribute nothing to the conditional
    # likelihood, so remove them.
    has_lower <- lower_total > 0

    if (!any(has_lower)) {
      stop(
        "No grouped/count rows have positive counts below `cutoff_level`; stage 1 cannot be fit.",
        call. = FALSE
      )
    }

    d1 <- cov_mf[has_lower, , drop = FALSE]

    y1_mat <- lower_counts[has_lower, , drop = FALSE]
    colnames(y1_mat) <- lower_levels

    d1$.tsco_y_stage1 <- I(y1_mat)

    # Stage 2: collapsed marginal model {Y < C, C, ..., K - 1}.
    d2 <- cov_mf

    y2_mat <- cbind(lower_total, upper_counts)
    colnames(y2_mat) <- collapsed_levels

    d2$.tsco_y_stage2 <- I(y2_mat)

    y1_counts <- y1_mat
    y2_counts <- y2_mat

    stage1_formula <- make_stage_formula(".tsco_y_stage1")
    stage2_formula <- make_stage_formula(".tsco_y_stage2")

    w1 <- if (is.null(w_all)) NULL else w_all[has_lower]
    w2 <- w_all

  } else {
    # Individual-level data.

    is_lower <- y_chr %in% lower_levels

    n_groups <- nrow(mf)
    n_obs <- nrow(mf)
    n_stage1_groups <- sum(is_lower)
    n_stage1_obs <- sum(is_lower)

    d1 <- cov_mf[is_lower, , drop = FALSE]

    if (stage1 == "po") {
      d1$.tsco_y_stage1 <- ordered(y_chr[is_lower], levels = lower_levels)
    } else {
      d1$.tsco_y_stage1 <- factor(y_chr[is_lower], levels = lower_levels)
    }

    d2 <- cov_mf

    y2_chr <- ifelse(y_chr %in% lower_levels, lower_collapsed_label, y_chr)

    if (stage2 == "po") {
      d2$.tsco_y_stage2 <- ordered(y2_chr, levels = collapsed_levels)
    } else {
      d2$.tsco_y_stage2 <- factor(y2_chr, levels = collapsed_levels)
    }

    y1_counts <- .tsco_factor_to_counts(
      d1$.tsco_y_stage1,
      levels = lower_levels
    )

    y2_counts <- .tsco_factor_to_counts(
      d2$.tsco_y_stage2,
      levels = collapsed_levels
    )

    stage1_formula <- make_stage_formula(".tsco_y_stage1")
    stage2_formula <- make_stage_formula(".tsco_y_stage2")

    w1 <- if (is.null(w_all)) NULL else w_all[is_lower]
    w2 <- w_all
  }

  # --------------------------------------------------------------------------
  # Fit models.
  # --------------------------------------------------------------------------

  fit1 <- fit_vglm(
    stage1_formula,
    make_family(stage1),
    d1,
    w1,
    ...
  )

  fit2 <- fit_vglm(
    stage2_formula,
    make_family(stage2),
    d2,
    w2,
    ...
  )

  p1_hat <- stats::predict(fit1, type = "response")
  p2_hat <- stats::predict(fit2, type = "response")

  p1_hat <- .tsco_align_prob(p1_hat, lower_levels)
  p2_hat <- .tsco_align_prob(p2_hat, collapsed_levels)

  if (!is.null(w1)) {
    y1_counts_ll <- y1_counts * as.numeric(w1)
  } else {
    y1_counts_ll <- y1_counts
  }

  if (!is.null(w2)) {
    y2_counts_ll <- y2_counts * as.numeric(w2)
  } else {
    y2_counts_ll <- y2_counts
  }

  logLik_stage1 <- .tsco_count_loglik(y1_counts_ll, p1_hat)
  logLik_stage2 <- .tsco_count_loglik(y2_counts_ll, p2_hat)
  logLik_total <- logLik_stage1 + logLik_stage2

  out <- list(
    call = match.call(),
    formula = formula,
    stage1_formula = stage1_formula,
    stage2_formula = stage2_formula,
    stage1 = stage1,
    stage2 = stage2,
    cutoff_level = cutoff_level,
    cutoff_index = cutoff_index,
    levels = y_levels,
    lower_levels = lower_levels,
    upper_levels = upper_levels,
    lower_collapsed_label = lower_collapsed_label,
    collapsed_levels = collapsed_levels,
    po.reverse = po.reverse,
    grouped = is_grouped,
    fit_stage1 = fit1,
    fit_stage2 = fit2,
    predict_data = d2,
    # sample-size accounting
    n = n_obs,
    n_obs = n_obs,
    n_groups = n_groups,
    n_stage1 = n_stage1_obs,
    n_stage1_obs = n_stage1_obs,
    n_stage1_groups = n_stage1_groups,
    logLik_stage1 = logLik_stage1,
    logLik_stage2 = logLik_stage2,
    logLik = logLik_total
  )


  class(out) <- "tsco"
  out
}
