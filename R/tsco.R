#' Fit a Two-Stage Conditional Odds (TSCO) Model
#'
#' Fits a two-stage conditional odds model for an ordinal outcome. Backend
#' selection is automatic: individual-level proportional-odds stages use
#' [rms::orm()], while all multinomial or grouped-count stages use
#' [VGAM::vglm()]. Set `po_engine = "orm"` or `po_engine = "vglm"` to
#' explicitly request an engine for proportional-odds stages.
#'
#' The outcome is partitioned at a user-specified cutoff C. Stage 1 models
#' Y | Y < C. Stage 2 models a collapsed outcome that combines all categories
#' below C with C, ..., K - 1.
#'
#' The response may be either an individual-level ordinal outcome or a grouped
#' count response supplied using cbind(), as in
#' cbind(normal, mild, severe) ~ x.
#'
#' @param formula Model formula, e.g. y ~ x1 + x2 for individual-level data or
#'   cbind(y0, y1, y2) ~ x1 + x2 for grouped count data. Inline
#'   transformations such as `log(x)` or `poly(x, 2)` are supported: the
#'   stage-specific right-hand sides are rebuilt from the model-frame columns,
#'   and `newdata` supplied to [predict.tsco()] may contain either the raw
#'   variables or the transformed columns.
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
#' @param stage2 Model for the collapsed marginal outcome formed by combining
#'   all levels below the cutoff with C, ..., K - 1. Either "multinomial" or
#'   "po".
#' @param po_engine Engine for proportional-odds stages. `"auto"` (the
#'   default) uses [rms::orm()] for individual-level responses and
#'   [VGAM::vglm()] for grouped-count responses. `"orm"` explicitly requests
#'   [rms::orm()] and errors for grouped-count PO stages; `"vglm"` always uses
#'   [VGAM::vglm()].
#' @param stage1_args Named list of additional arguments passed only to the
#'   selected stage-1 fitter. Names must be unique and cannot be `formula`,
#'   `data`, `family`, `weights`, or `na.action`, which are controlled by
#'   `tsco()`. These arguments are reused for reduced-model likelihood-ratio
#'   test refits, so arguments tied to the full model's parameter dimension may
#'   prevent those tests from being calculated.
#' @param stage2_args Named list of additional arguments passed only to the
#'   selected stage-2 fitter, with the same restrictions and likelihood-ratio
#'   refit behavior as `stage1_args`.
#' @param weights Optional finite, nonnegative multiplicative case weights with
#'   at least one positive value. For grouped responses, a weight multiplies the
#'   complete count contribution from that covariate pattern.
#' @param na.action Missing-data action passed to model.frame().
#' @param warn_degenerate If TRUE, warn when K = 3, which reduces to two binary
#'   regressions.
#'
#' @section Limitations:
#' Prediction-time stage-1 support checks compare factor levels marginally and
#' do not detect absent combinations of observed levels induced by interaction
#' terms. See [predict.tsco()] for details.
#'
#' Complete or quasi-complete separation in either fitted stage can yield
#' infinite or unstable maximum-likelihood estimates. Likelihood-ratio tests do
#' not remove this problem, and their ordinary chi-squared reference
#' distribution may be unreliable under separation. [summary.tsco()] reports a
#' diagnostic when a stage shows the usual symptoms, but the package does not
#' correct separation.
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
    po_engine = c("auto", "orm", "vglm"),
    stage1_args = list(),
    stage2_args = list(),
    weights = NULL,
    na.action = stats::na.omit,
    warn_degenerate = TRUE) {

  stage1 <- match.arg(stage1)
  stage2 <- match.arg(stage2)
  po_engine <- match.arg(po_engine)
  stage1_args <- .tsco_validate_stage_args(stage1_args, "stage1_args")
  stage2_args <- .tsco_validate_stage_args(stage2_args, "stage2_args")

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

  # Convert character predictors to factors using the full analysis dataset.
  # This prevents character predictors from being re-leveled differently in
  # stage 1 and stage 2 after subsetting.
  char_predictors <- names(cov_mf)[
    vapply(cov_mf, is.character, logical(1L))
  ]

  for (nm in char_predictors) {
    cov_mf[[nm]] <- factor(
      cov_mf[[nm]],
      levels = sort(unique(cov_mf[[nm]]))
    )
  }
  predict_data <- cov_mf

  # Store factor metadata for prediction-time newdata preparation.
  factor_predictors <- names(cov_mf)[
    vapply(cov_mf, is.factor, logical(1L))
  ]

  predictor_levels <- lapply(cov_mf[factor_predictors], levels)
  predictor_ordered <- lapply(cov_mf[factor_predictors], is.ordered)

  # Levels actually observed in the full analysis data. Factor levels that are
  # declared but never observed cannot support prediction.
  predictor_observed_levels <- lapply(
    cov_mf[factor_predictors],
    function(z) {
      lev <- levels(z)
      lev[lev %in% as.character(z)]
    }
  )

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
    if (!is.numeric(weights) || !is.null(dim(weights))) {
      stop("`weights` must be a numeric vector.", call. = FALSE)
    }

    if (length(weights) == 0L ||
        any(!is.finite(weights)) ||
        any(weights < 0) ||
        !any(weights > 0)) {
      stop(
        "`weights` must be finite, nonnegative, and contain at least one positive value.",
        call. = FALSE
      )
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

  # Internal helper: construct a stage formula with a new response.
  #
  # The right-hand side is rebuilt from the model-frame column names rather
  # than reused verbatim, because the stage data frames below hold evaluated
  # model-frame columns. For a formula of plain variables this reproduces the
  # original right-hand side exactly; for one containing inline
  # transformations it resolves them against the columns that exist.
  rhs_labels <- .tsco_stage_rhs_labels(tt)

  make_stage_formula <- function(response_name) {
    rhs <- if (length(rhs_labels) == 0L) {
      "1"
    } else {
      paste(rhs_labels, collapse = " + ")
    }

    if (identical(as.integer(attr(tt, "intercept")), 0L)) {
      rhs <- paste(rhs, "- 1")
    }

    stats::as.formula(
      paste0("`", response_name, "` ~ ", rhs),
      env = environment(formula)
    )
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

  # Remove stage-1 factor levels with no contributing observations. Retaining
  # their all-zero design columns can make orm fits singular. Full-data factor
  # metadata remain stored separately for prediction support checks.
  for (nm in intersect(factor_predictors, names(d1))) {
    d1[[nm]] <- droplevels(d1[[nm]])
  }

  if (!is.null(w1) && !any(w1 > 0)) {
    stop(
      "`weights` must contain at least one positive value among observations contributing to stage 1.",
      call. = FALSE
    )
  }

  # Factor levels actually represented in the stage-1 fitting data.
  # If a factor level is absent from stage 1, conditional probabilities
  # Y | Y < C for that level are not identified from the stage-1 likelihood.
  stage1_factor_levels <- lapply(
    factor_predictors,
    function(nm) {
      if (nm %in% names(d1)) {
        z <- d1[[nm]]

        if (is.factor(z)) {
          lev <- levels(z)
          lev[lev %in% as.character(z)]
        } else {
          sort(unique(as.character(z)))
        }
      } else {
        character()
      }
    }
  )

  names(stage1_factor_levels) <- factor_predictors

  # --------------------------------------------------------------------------
  # Fit models. `orm()` is used only where it supports the response structure:
  # individual-level proportional-odds stages. All other stages use VGAM.
  # --------------------------------------------------------------------------

  stage1_engine <- .tsco_stage_engine(
    stage1,
    grouped = is_grouped,
    po_engine = po_engine,
    weighted = !is.null(w1)
  )
  stage2_engine <- .tsco_stage_engine(
    stage2,
    grouped = is_grouped,
    po_engine = po_engine,
    weighted = !is.null(w2)
  )

  fit1 <- .tsco_fit_stage(
    formula = stage1_formula,
    kind = stage1,
    engine = stage1_engine,
    data = d1,
    weights = w1,
    stage_args = stage1_args
  )

  fit2 <- .tsco_fit_stage(
    formula = stage2_formula,
    kind = stage2,
    engine = stage2_engine,
    data = d2,
    weights = w2,
    stage_args = stage2_args
  )

  p1_hat <- .tsco_predict_stage_prob(
    fit1,
    engine = stage1_engine,
    newdata = d1,
    levels = lower_levels
  )
  p2_hat <- .tsco_predict_stage_prob(
    fit2,
    engine = stage2_engine,
    newdata = d2,
    levels = collapsed_levels
  )

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
    grouped = is_grouped,
    po_engine = po_engine,
    stage1_engine = stage1_engine,
    stage2_engine = stage2_engine,
    fit_stage1 = fit1,
    fit_stage2 = fit2,

    # For prediction and LRT refits
    terms_rhs = stats::delete.response(tt),
    predict_data = predict_data,
    data_stage1 = d1,
    data_stage2 = d2,
    weights_stage1 = w1,
    weights_stage2 = w2,

    # Weighted category totals, used for closed-form intercept-only fits in
    # joint likelihood-ratio tests.
    totals_stage1 = colSums(y1_counts_ll),
    totals_stage2 = colSums(y2_counts_ll),
    stage1_args = stage1_args,
    stage2_args = stage2_args,

    # Predictor metadata for prediction
    predictor_levels = predictor_levels,
    predictor_ordered = predictor_ordered,
    predictor_observed_levels = predictor_observed_levels,
    stage1_factor_levels = stage1_factor_levels,

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
