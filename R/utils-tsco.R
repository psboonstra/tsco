# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.tsco_deparse <- function(x, width.cutoff = 500L) {
  paste(deparse(x, width.cutoff = width.cutoff), collapse = " ")
}

.tsco_df <- function(fit) {
  valid_df <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  }

  if (inherits(fit, "orm")) {
    if (isTRUE(fit$fail)) {
      return(NA_real_)
    }

    df <- tryCatch(
      attr(stats::logLik(fit), "df", exact = TRUE),
      error = function(e) NULL
    )

    if (valid_df(df)) {
      return(as.numeric(df))
    }
  }

  if (inherits(fit, "vglm")) {
    rank <- tryCatch(
      methods::slot(fit, "rank"),
      error = function(e) NULL
    )

    if (valid_df(rank)) {
      return(as.numeric(rank))
    }
  }

  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)

  if (is.null(cf)) {
    return(NA_real_)
  }

  df <- tryCatch(sum(is.finite(cf)), error = function(e) NA_real_)
  if (valid_df(df)) as.numeric(df) else NA_real_
}

.tsco_stage_vcov <- function(fit) {
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)

  if (is.null(cf)) {
    return(NULL)
  }

  V <- if (inherits(fit, "orm")) {
    tryCatch(
      stats::vcov(fit, intercepts = "all"),
      error = function(e) {
        tryCatch(stats::vcov(fit), error = function(e) NULL)
      }
    )
  } else {
    tryCatch(stats::vcov(fit), error = function(e) NULL)
  }

  if (is.null(V)) {
    return(NULL)
  }

  V <- as.matrix(V)

  if (is.null(rownames(V)) || is.null(colnames(V))) {
    if (nrow(V) == length(cf) && ncol(V) == length(cf)) {
      dimnames(V) <- list(names(cf), names(cf))
    }

    return(V)
  }

  if (identical(rownames(V), names(cf)) && identical(colnames(V), names(cf))) {
    return(V)
  }

  # A reduced orm covariance may contain only the middle threshold and slopes.
  # Preserve the full layout while marking unavailable covariances as NA.
  out <- matrix(
    NA_real_,
    nrow = length(cf),
    ncol = length(cf),
    dimnames = list(names(cf), names(cf))
  )
  common <- intersect(names(cf), rownames(V))
  out[common, common] <- V[common, common, drop = FALSE]
  out
}


.tsco_stage_coef_table <- function(
    fit,
    kind,
    engine = NULL) {
  b_raw <- tryCatch(stats::coef(fit), error = function(e) NULL)
  V_raw <- .tsco_stage_vcov(fit)

  if (is.null(b_raw)) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  s <- .tsco_sign_vector(
    fit,
    kind = kind,
    engine = engine
  )

  b <- b_raw * s

  if (is.null(V_raw)) {
    return(
      matrix(
        b,
        ncol = 1L,
        dimnames = list(names(b), "Estimate")
      )
    )
  }

  V_raw <- as.matrix(V_raw)

  if (!is.null(rownames(V_raw))) {
    s_v <- .tsco_align_sign_vector(s, V_raw)
    V <- V_raw * outer(s_v, s_v)
    b <- b[rownames(V_raw)]
  } else {
    V <- V_raw * outer(s, s)
  }

  se <- sqrt(diag(V))

  z <- b / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

  out <- cbind(
    Estimate = b,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = p
  )

  rownames(out) <- names(b)
  out
}

#' Construct a stage-level coefficient table on the reported TSCO scale.
#'
#' Estimates and covariance matrices are transformed using `.tsco_sign_vector()`.
#' Standard errors and p-values are computed from the transformed covariance
#' matrix.
#'
#' @param fit A fitted stage model (a "vglm" object).
#' @param kind "po" or "multinomial", i.e. object$stage1 or object$stage2.
#' @param engine The stage fitting backend, either "orm" or "vglm".
.tsco_coef_table <- function(fit, kind, engine = NULL) {
  .tsco_stage_coef_table(
    fit = fit,
    kind = kind,
    engine = engine
  )
}

.tsco_model_label <- function(x) {
  switch(
    x,
    po = "proportional odds",
    multinomial = "multinomial",
    x
  )
}

.tsco_count_loglik <- function(y, p, eps = .Machine$double.eps) {
  y <- as.matrix(y)
  p <- as.matrix(p)

  if (nrow(y) != nrow(p)) {
    stop(
      "Internal error: response and probability matrices have different numbers of rows.",
      call. = FALSE
    )
  }

  if (ncol(y) != ncol(p)) {
    stop(
      "Internal error: response and probability matrices have different numbers of columns.",
      call. = FALSE
    )
  }

  if (any(!is.finite(y))) {
    stop(
      "Internal error: response count matrix contains non-finite values.",
      call. = FALSE
    )
  }

  if (any(y < 0)) {
    stop(
      "Internal error: response count matrix contains negative values.",
      call. = FALSE
    )
  }

  if (any(!is.finite(p))) {
    stop(
      "Internal error: predicted probability matrix contains non-finite values.",
      call. = FALSE
    )
  }

  if (any(p < 0)) {
    stop(
      "Internal error: predicted probability matrix contains negative values.",
      call. = FALSE
    )
  }

  p <- pmax(p, eps)
  p <- p / rowSums(p)

  sum(y * log(p))
}

.tsco_factor_to_counts <- function(y, levels) {
  y_chr <- as.character(y)

  idx <- match(y_chr, levels)

  if (anyNA(idx)) {
    bad <- unique(y_chr[is.na(idx)])
    stop(
      "Internal error: outcome values not found in supplied levels: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  out <- matrix(
    0,
    nrow = length(y_chr),
    ncol = length(levels),
    dimnames = list(NULL, levels)
  )

  out[cbind(seq_along(y_chr), idx)] <- 1
  out
}

#' Sign vector correcting a fitted stage's coefficients to the reported TSCO
#' scale.
#'
#' The manuscript parameterizes both submodels with a subtractive convention
#' on x'beta (Eq. 1 for PO, Eq. 2 for MR/multinomial). Multinomial coefficients
#' are reported in that parameterization. For PO stages, covariate coefficients
#' use the manuscript's sign convention, but thresholds use the common
#' lower-tail representation `Pr(Y <= j | X) = expit(gamma_j + x'beta)`.
#' Thus `gamma_j = -alpha_(j + 1)` for manuscript upper-tail cutpoints. The
#' required transformation depends on the selected fitting backend:
#'   - "po" with `VGAM::cumulative(reverse = FALSE)` already uses the reported
#'     lower-tail representation, so no coefficients are flipped.
#'   - "po" with `rms::orm()` uses the opposite-sign upper-tail representation.
#'     Negating the entire coefficient vector converts both thresholds and
#'     slopes to the reported lower-tail representation.
#'   - "multinomial" with `VGAM::multinomial(refLevel = 1)` uses the additive
#'     convention `eta = alpha + x'beta`. Every non-intercept coefficient is
#'     negated to obtain the manuscript's `alpha - x'beta` convention.
#'
#' This helper is the single point of control for coefficient transformations,
#' keeping the sign logic consistent across `coef.tsco()`, `vcov.tsco()`, and
#' `summary.tsco()`.
#'
#' @param fit A fitted stage model.
#' @param kind "po" or "multinomial", i.e. object$stage1 or object$stage2.
#' @return A named numeric vector of +1/-1, same length and names as
#'   stats::coef(fit).
#' @param engine The stage fitting backend, either "orm" or "vglm".
.tsco_sign_vector <- function(
    fit,
    kind,
    engine = NULL) {

  cf <- stats::coef(fit)
  s <- rep(1, length(cf))
  names(s) <- names(cf)

  if (is.null(engine)) {
    engine <- if (inherits(fit, "orm")) "orm" else "vglm"
  }

  is_intercept <- if (identical(engine, "orm")) {
    seq_along(cf) <= fit$non.slopes
  } else {
    grepl("(Intercept)", names(cf), fixed = TRUE)
  }

  if (identical(kind, "multinomial")) {
    # VGAM multinomial uses alpha + x'beta while the manuscript uses
    # alpha - x'beta.
    s[!is_intercept] <- -1
  }

  if (identical(engine, "orm")) {
    # Convert orm's upper-tail coefficients to the common lower-tail scale.
    s[] <- -1
  }

  s
}

.tsco_align_sign_vector <- function(s, V) {
  rn <- rownames(V)

  if (is.null(rn)) {
    return(s)
  }

  if (!all(rn %in% names(s))) {
    stop(
      "Internal error: covariance matrix names do not align with coefficient names.",
      call. = FALSE
    )
  }

  s[rn]
}

.tsco_align_prob <- function(p, levels) {
  if (is.null(dim(p))) {
    p <- matrix(p, nrow = 1L)
  } else {
    p <- as.matrix(p)
  }

  if (ncol(p) != length(levels)) {
    stop(
      "Predicted probability matrix has ", ncol(p),
      " columns, but ", length(levels), " levels were expected.",
      call. = FALSE
    )
  }

  cn <- colnames(p)

  if (is.null(cn) || anyNA(cn) || any(cn == "")) {
    colnames(p) <- levels
    return(p[, levels, drop = FALSE])
  }

  missing <- setdiff(levels, cn)

  if (length(missing) == 0L) {
    return(p[, levels, drop = FALSE])
  }

  # Last-resort positional fallback. This handles cases where VGAM returns
  # probabilities in the correct order but with transformed/syntactic names.
  colnames(p) <- levels
  p[, levels, drop = FALSE]
}


.tsco_base_coef_name <- function(x) {
  # Remove package-level stage prefix from coef.tsco().
  x <- sub("^stage[12]:", "", x)

  # Collapse RMS threshold names so they are excluded from joint tests.
  x[x == "Intercept" | grepl(">=", x, fixed = TRUE)] <- "(Intercept)"

  # Remove VGAM linear-predictor suffix, e.g. let:1 -> let.
  # This preserves interactions such as x:z while converting x:z:1 to x:z.
  x <- sub(":[0-9]+$", "", x)

  x
}

.tsco_is_intercept_name <- function(x) {
  x == "(Intercept)"
}

.tsco_wald_stat <- function(b, V) {
  b <- as.numeric(b)
  V <- as.matrix(V)

  if (length(b) == 0L) {
    return(NA_real_)
  }

  if (any(!is.finite(b)) || any(!is.finite(V))) {
    return(NA_real_)
  }

  tryCatch(
    as.numeric(crossprod(b, qr.solve(V, b))),
    error = function(e) NA_real_
  )
}

.tsco_stage_coef_terms <- function(fit, formula, data, engine, stage) {
  cf <- stats::coef(fit)
  coef_names <- names(cf)

  if (is.null(coef_names)) {
    stop(
      "Cannot map coefficients to formula terms because coefficients are unnamed.",
      call. = FALSE
    )
  }

  tt <- stats::terms(formula)
  term_labels <- attr(tt, "term.labels")
  mm <- stats::model.matrix(stats::delete.response(tt), data = data)
  assign <- attr(mm, "assign")
  design_idx <- which(assign > 0L)
  design_names <- colnames(mm)[design_idx]
  design_terms <- term_labels[assign[design_idx]]

  out <- rep(NA_character_, length(cf))
  names(out) <- paste0(stage, ":", coef_names)

  if (identical(engine, "orm")) {
    slope_idx <- if (fit$non.slopes < length(cf)) {
      seq.int(fit$non.slopes + 1L, length(cf))
    } else {
      integer()
    }

    if (length(slope_idx) != length(design_terms)) {
      stop(
        "Cannot align orm coefficients with formula terms.",
        call. = FALSE
      )
    }

    out[slope_idx] <- design_terms
    return(out)
  }

  base_names <- .tsco_base_coef_name(coef_names)
  slope_idx <- which(!.tsco_is_intercept_name(base_names))
  matched <- match(base_names[slope_idx], design_names)

  if (anyNA(matched)) {
    stop(
      "Cannot align vglm coefficients with formula terms: ",
      paste(unique(base_names[slope_idx][is.na(matched)]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  out[slope_idx] <- design_terms[matched]
  out
}

.tsco_joint_wald_tests <- function(object) {
  b <- stats::coef(object, stage = "both")
  V <- stats::vcov(object, stage = "both")

  if (length(b) == 0L) {
    return(data.frame())
  }

  V <- as.matrix(V)

  if (is.null(names(b))) {
    stop(
      "Cannot compute joint tests because coefficients are unnamed.",
      call. = FALSE
    )
  }

  if (is.null(rownames(V)) || is.null(colnames(V))) {
    stop(
      "Cannot compute joint tests because covariance matrix is unnamed.",
      call. = FALSE
    )
  }

  if (!all(names(b) %in% rownames(V)) || !all(names(b) %in% colnames(V))) {
    stop(
      "Cannot compute joint tests because coefficient and covariance names do not align.",
      call. = FALSE
    )
  }

  coef_terms <- c(
    .tsco_stage_coef_terms(
      object$fit_stage1,
      object$stage1_formula,
      object$data_stage1,
      object$stage1_engine,
      "stage1"
    ),
    .tsco_stage_coef_terms(
      object$fit_stage2,
      object$stage2_formula,
      object$data_stage2,
      object$stage2_engine,
      "stage2"
    )
  )

  if (!all(names(b) %in% names(coef_terms))) {
    stop(
      "Cannot compute joint tests because coefficients could not be mapped to formula terms.",
      call. = FALSE
    )
  }

  coef_terms <- coef_terms[names(b)]
  term_labels <- attr(stats::terms(object$stage1_formula), "term.labels")
  test_names <- term_labels[term_labels %in% unique(coef_terms[!is.na(coef_terms)])]

  if (length(test_names) == 0L) {
    return(data.frame())
  }

  res <- lapply(test_names, function(term) {
    idx <- which(coef_terms == term)

    b_term <- b[idx]
    V_term <- V[names(b_term), names(b_term), drop = FALSE]

    stat <- .tsco_wald_stat(b_term, V_term)
    df <- length(b_term)

    p <- if (is.finite(stat)) {
      stats::pchisq(stat, df = df, lower.tail = FALSE)
    } else {
      NA_real_
    }

    data.frame(
      term = term,
      df = df,
      Chisq = stat,
      `Pr(>Chisq)` = p,
      check.names = FALSE
    )
  })

  out <- do.call(rbind, res)
  rownames(out) <- out$term
  out$term <- NULL

  out
}

.tsco_stage_engine <- function(kind, grouped, po_engine) {
  if (!identical(kind, "po")) {
    return("vglm")
  }

  if (isTRUE(grouped)) {
    if (identical(po_engine, "orm")) {
      stop(
        "`po_engine = \"orm\"` cannot be used with grouped-count PO stages; ",
        "use `po_engine = \"auto\"` or `po_engine = \"vglm\"`.",
        call. = FALSE
      )
    }

    return("vglm")
  }

  if (identical(po_engine, "vglm")) "vglm" else "orm"
}

.tsco_validate_stage_args <- function(args, arg_name) {
  if (!is.list(args)) {
    stop("`", arg_name, "` must be a list.", call. = FALSE)
  }

  if (length(args) == 0L) {
    return(args)
  }

  nms <- names(args)

  if (is.null(nms) || anyNA(nms) || any(nms == "")) {
    stop("Every element of `", arg_name, "` must be named.", call. = FALSE)
  }

  duplicated_names <- unique(nms[duplicated(nms)])

  if (length(duplicated_names) > 0L) {
    stop(
      "Duplicated argument names in `", arg_name, "`: ",
      paste(duplicated_names, collapse = ", "), ".",
      call. = FALSE
    )
  }

  reserved <- intersect(
    c("formula", "data", "family", "weights", "na.action"),
    nms
  )

  if (length(reserved) > 0L) {
    stop(
      "`", arg_name, "` cannot contain arguments controlled by `tsco()`: ",
      paste(reserved, collapse = ", "), ".",
      call. = FALSE
    )
  }

  args
}

.tsco_fit_stage <- function(
    formula,
    kind,
    engine,
    data,
    weights = NULL,
    stage_args = list()) {

  if (identical(engine, "orm")) {
    if (!identical(kind, "po")) {
      stop("The orm backend supports only proportional-odds stages.", call. = FALSE)
    }

    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("Package 'rms' is required for individual-level PO stages.", call. = FALSE)
    }

    args <- list(
      formula = formula,
      data = data,
      family = "logistic"
    )

    if (!is.null(weights)) {
      args$weights <- weights
    }

    return(do.call(rms::orm, c(args, stage_args)))
  }

  if (!identical(engine, "vglm")) {
    stop("Unknown stage fitting engine: ", engine, call. = FALSE)
  }

  if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop("Package 'VGAM' is required for this stage.", call. = FALSE)
  }

  family <- if (identical(kind, "po")) {
    VGAM::cumulative(
      link = "logitlink",
      parallel = TRUE,
      reverse = FALSE
    )
  } else if (identical(kind, "multinomial")) {
    VGAM::multinomial(refLevel = 1)
  } else {
    stop("Unknown stage model kind: ", kind, call. = FALSE)
  }

  args <- list(
    formula = formula,
    family = family,
    data = data
  )

  if (!is.null(weights)) {
    args$weights <- weights
  }

  do.call(VGAM::vglm, c(args, stage_args))
}

.tsco_predict_stage_prob <- function(
    fit,
    engine,
    newdata,
    levels = NULL) {

  if (identical(engine, "orm")) {
    out <- stats::predict(
      fit,
      newdata = newdata,
      type = "fitted.ind"
    )

    if (is.null(dim(out)) && !is.null(levels) && length(levels) == 2L) {
      out <- cbind(1 - out, out)
      colnames(out) <- levels
    }

    return(out)
  }

  if (identical(engine, "vglm")) {
    return(VGAM::predictvglm(
      fit,
      newdata = newdata,
      type = "response"
    ))
  }

  stop("Unknown stage fitting engine: ", engine, call. = FALSE)
}

.tsco_stage_loglik <- function(fit, engine) {
  if (identical(engine, "orm")) {
    out <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
  } else if (identical(engine, "vglm")) {
    out <- tryCatch({
      as.numeric(fit@criterion["loglikelihood"])
    }, error = function(e) NA_real_)
  } else {
    out <- NA_real_
  }

  if (length(out) != 1L || !is.finite(out)) NA_real_ else out
}

.tsco_drop_term_formula <- function(formula, term) {
  tt <- stats::terms(formula)
  term_labels <- attr(tt, "term.labels")

  pos <- match(term, term_labels)

  if (is.na(pos)) {
    stop("Term not found in formula: ", term, call. = FALSE)
  }

  tt_reduced <- stats::drop.terms(
    tt,
    dropx = pos,
    keep.response = TRUE
  )

  f_reduced <- stats::formula(tt_reduced)
  environment(f_reduced) <- environment(formula)

  f_reduced
}

.tsco_clean_lrt <- function(lrt, term, tolerance = 1e-7) {
  if (!is.finite(lrt)) {
    return(list(
      value = NA_real_,
      problem = paste0(
        "Term `", term, "`: the likelihood-ratio statistic was not finite."
      )
    ))
  }

  if (lrt <= -tolerance) {
    return(list(
      value = NA_real_,
      problem = paste0(
        "Term `", term, "`: the reduced model had a larger log-likelihood ",
        "than the full model; the likelihood-ratio statistic was set to NA."
      )
    ))
  }

  if (lrt < 0) {
    lrt <- 0
  }

  list(value = lrt, problem = NA_character_)
}

.tsco_joint_lrt_tests <- function(object) {
  if (is.null(object$data_stage1) ||
      is.null(object$data_stage2) ||
      is.null(object$weights_stage1) && !("weights_stage1" %in% names(object)) ||
      is.null(object$weights_stage2) && !("weights_stage2" %in% names(object))) {
    stop(
      "Joint LRTs require stage-specific data stored in the tsco object.",
      call. = FALSE
    )
  }

  term_labels <- attr(stats::terms(object$stage1_formula), "term.labels")

  if (length(term_labels) == 0L) {
    return(
      data.frame(
        df = integer(),
        Chisq = numeric(),
        `Pr(>Chisq)` = numeric(),
        check.names = FALSE
      )
    )
  }

  stage1_engine <- object$stage1_engine
  stage2_engine <- object$stage2_engine

  ll1_full <- .tsco_stage_loglik(object$fit_stage1, stage1_engine)
  ll2_full <- .tsco_stage_loglik(object$fit_stage2, stage2_engine)

  df1_full <- .tsco_df(object$fit_stage1)
  df2_full <- .tsco_df(object$fit_stage2)

  if (!is.finite(ll1_full) || !is.finite(ll2_full)) {
    stop(
      "Could not extract full-model log-likelihoods for LRTs.",
      call. = FALSE
    )
  }

  stage1_args <- object$stage1_args
  stage2_args <- object$stage2_args

  res <- lapply(term_labels, function(term) {
    tryCatch({
      f1_reduced <- .tsco_drop_term_formula(object$stage1_formula, term)
      f2_reduced <- .tsco_drop_term_formula(object$stage2_formula, term)

      fit1_reduced <- .tsco_fit_stage(
        formula = f1_reduced,
        kind = object$stage1,
        engine = stage1_engine,
        data = object$data_stage1,
        weights = object$weights_stage1,
        stage_args = stage1_args
      )

      fit2_reduced <- .tsco_fit_stage(
        formula = f2_reduced,
        kind = object$stage2,
        engine = stage2_engine,
        data = object$data_stage2,
        weights = object$weights_stage2,
        stage_args = stage2_args
      )

      ll1_reduced <- .tsco_stage_loglik(fit1_reduced, stage1_engine)
      ll2_reduced <- .tsco_stage_loglik(fit2_reduced, stage2_engine)

      if (!is.finite(ll1_reduced) || !is.finite(ll2_reduced)) {
        stop("Could not extract reduced-model log-likelihoods.", call. = FALSE)
      }

      df1_reduced <- .tsco_df(fit1_reduced)
      df2_reduced <- .tsco_df(fit2_reduced)

      lrt <- 2 * (
        (ll1_full + ll2_full) -
          (ll1_reduced + ll2_reduced)
      )
      cleaned <- .tsco_clean_lrt(lrt, term)
      lrt <- cleaned$value

      df <- (df1_full - df1_reduced) + (df2_full - df2_reduced)

      p <- if (is.finite(lrt) && is.finite(df) && df > 0) {
        stats::pchisq(lrt, df = df, lower.tail = FALSE)
      } else {
        NA_real_
      }

      list(
        row = data.frame(
          term = term,
          df = df,
          Chisq = lrt,
          `Pr(>Chisq)` = p,
          check.names = FALSE
        ),
        error = NA_character_,
        problem = cleaned$problem
      )
    }, error = function(e) {
      message <- conditionMessage(e)

      list(
        row = data.frame(
          term = term,
          df = NA_real_,
          Chisq = NA_real_,
          `Pr(>Chisq)` = NA_real_,
          check.names = FALSE
        ),
        error = message,
        problem = paste0(
          "Term `", term, "`: reduced-model fitting failed: ", message
        )
      )
    })
  })

  problems <- vapply(res, function(x) x$problem, character(1L))
  problems <- problems[!is.na(problems)]

  if (length(problems) > 0L) {
    warning(
      paste(
        "Joint likelihood-ratio test problems:",
        paste0("- ", problems, collapse = "\n"),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  out <- do.call(rbind, lapply(res, function(x) x$row))
  rownames(out) <- out$term
  out$term <- NULL

  errors <- vapply(res, function(x) x$error, character(1L))
  names(errors) <- term_labels
  errors <- errors[!is.na(errors)]

  if (length(errors) > 0L) {
    attr(out, "errors") <- errors
  }

  out
}

.tsco_prepare_newdata <- function(object, newdata) {
  if (is.null(newdata)) {
    if (is.null(object$predict_data)) {
      stop(
        "No prediction data stored in `object`; please supply `newdata`.",
        call. = FALSE
      )
    }

    return(object$predict_data)
  }

  nd <- as.data.frame(newdata)

  predictor_levels <- object$predictor_levels
  predictor_ordered <- object$predictor_ordered
  predictor_observed_levels <- object$predictor_observed_levels

  if (is.null(predictor_levels) || length(predictor_levels) == 0L) {
    return(nd)
  }

  for (nm in names(predictor_levels)) {
    if (!nm %in% names(nd)) {
      stop(
        "Prediction data are missing predictor `", nm, "`.",
        call. = FALSE
      )
    }

    lev <- predictor_levels[[nm]]
    ord <- isTRUE(predictor_ordered[[nm]])

    observed <- if (!is.null(predictor_observed_levels) &&
                    nm %in% names(predictor_observed_levels)) {
      predictor_observed_levels[[nm]]
    } else {
      lev
    }

    vals <- as.character(nd[[nm]])

    # Values never observed in the fitting data should not be predicted.
    bad <- setdiff(unique(vals[!is.na(vals)]), observed)

    if (length(bad) > 0L) {
      stop(
        "Prediction data contain level(s) for predictor `", nm,
        "` that were not observed in the fitting data: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }

    # Reconstruct factor using the original training-data level order.
    nd[[nm]] <- factor(vals, levels = lev, ordered = ord)
  }

  nd
}

.tsco_stage1_unsupported_rows <- function(object, newdata) {
  if (is.null(object$stage1_factor_levels) ||
      length(object$stage1_factor_levels) == 0L) {
    return(rep(FALSE, nrow(newdata)))
  }

  unsupported <- rep(FALSE, nrow(newdata))

  for (nm in names(object$stage1_factor_levels)) {
    if (!nm %in% names(newdata)) {
      next
    }

    observed_stage1 <- object$stage1_factor_levels[[nm]]
    vals <- as.character(newdata[[nm]])

    # Missing values or levels absent from the stage-1 fitting data imply that
    # conditional lower-partition probabilities are not identified.
    unsupported_nm <- is.na(vals) | !(vals %in% observed_stage1)

    unsupported <- unsupported | unsupported_nm
  }

  unsupported
}


.tsco_empirical_lower_probs <- function(object) {
  if (is.null(object$data_stage1) || is.null(object$stage1_formula)) {
    stop(
      "Cannot compute empirical lower-category probabilities because ",
      "stage-1 data are not stored in the TSCO object.",
      call. = FALSE
    )
  }

  lhs_vars <- all.vars(object$stage1_formula[[2L]])

  if (length(lhs_vars) != 1L) {
    stop(
      "Could not identify the stage-1 response variable.",
      call. = FALSE
    )
  }

  yname <- lhs_vars[1L]

  if (!yname %in% names(object$data_stage1)) {
    stop(
      "Stage-1 response variable not found in stored stage-1 data.",
      call. = FALSE
    )
  }

  y <- object$data_stage1[[yname]]

  if (is.matrix(y) || is.data.frame(y)) {
    counts <- colSums(as.matrix(y))
    counts <- counts[object$lower_levels]
  } else {
    counts <- table(
      factor(
        as.character(y),
        levels = object$lower_levels
      )
    )
    counts <- as.numeric(counts)
    names(counts) <- object$lower_levels
  }

  if (anyNA(counts) || sum(counts) <= 0) {
    warning(
      "Could not compute empirical lower-category probabilities; ",
      "using a uniform split instead.",
      call. = FALSE
    )

    out <- rep(
      1 / length(object$lower_levels),
      length(object$lower_levels)
    )
    names(out) <- object$lower_levels
    return(out)
  }

  out <- as.numeric(counts) / sum(counts)
  names(out) <- object$lower_levels

  out
}
