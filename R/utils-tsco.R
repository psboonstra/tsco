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
    engine = NULL,
    reported_names = NULL) {
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
    # Same reported names as the full table, `coef()` and `vcov()`; a fit
    # whose covariance is unavailable must not fall back to backend names.
    return(
      matrix(
        b,
        ncol = 1L,
        dimnames = list(
          .tsco_apply_reported_names(names(b), reported_names),
          "Estimate"
        )
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

  rownames(out) <- .tsco_apply_reported_names(names(b), reported_names)
  out
}

#' Look up reported names for a (possibly reordered) subset of coefficients.
#'
#' @param coef_names The backend's coefficient names, in the order wanted.
#' @param reported_names A character vector of reported names, named by the
#'   backend's coefficient names, as returned by
#'   `.tsco_stage_reported_names()`.
#' @noRd
.tsco_apply_reported_names <- function(coef_names, reported_names) {
  if (is.null(reported_names) || is.null(names(reported_names))) {
    return(coef_names)
  }

  out <- reported_names[coef_names]

  if (anyNA(out)) {
    stop(
      "Internal error: could not map coefficient names to reported names: ",
      paste(coef_names[is.na(out)], collapse = ", "),
      call. = FALSE
    )
  }

  unname(out)
}

#' Construct a stage-level coefficient table on the reported TSCO scale.
#'
#' Estimates and covariance matrices are transformed using `.tsco_sign_vector()`.
#' Standard errors and p-values are computed from the transformed covariance
#' matrix.
#'
#' @param fit A fitted stage model (a "vglm" or "orm" object).
#' @param kind "po" or "multinomial", i.e. object$stage1 or object$stage2.
#' @param engine The stage fitting backend, either "orm" or "vglm".
#' @param reported_names Optional character vector of reported coefficient
#'   names, named by the backend's own coefficient names, as returned by
#'   `.tsco_stage_reported_names()`.
#' @noRd
.tsco_coef_table <- function(
    fit,
    kind,
    engine = NULL,
    reported_names = NULL) {

  .tsco_stage_coef_table(
    fit = fit,
    kind = kind,
    engine = engine,
    reported_names = reported_names
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
#' on x'beta, for both the PO and the MR/multinomial submodel. Multinomial coefficients
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
#' @noRd
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

  # `rms::predict.orm(type = "fitted.ind")` prefixes every column with the
  # response name, as in `y=a`. Strip a common `<name>=` prefix so these match
  # by name rather than falling through to column order.
  stripped <- sub("^[^=]*=", "", cn)

  if (!anyDuplicated(stripped) && all(levels %in% stripped)) {
    p <- p[, match(levels, stripped), drop = FALSE]
    colnames(p) <- levels
    return(p)
  }

  # A backend may return the columns in the right order but with syntactically
  # transformed names. Match on syntactic names before giving up on names
  # altogether.
  safe_levels <- make.names(levels)
  safe_cn <- make.names(cn)

  if (!anyDuplicated(safe_levels) &&
      !anyDuplicated(safe_cn) &&
      all(safe_levels %in% safe_cn)) {

    p <- p[, match(safe_levels, safe_cn), drop = FALSE]
    colnames(p) <- levels
    return(p)
  }

  # Last-resort positional fallback. Silently trusting column order would
  # mis-assign probabilities without any error, so warn.
  warning(
    "Could not match predicted probability column names to outcome levels (",
    paste(cn, collapse = ", "), " vs ", paste(levels, collapse = ", "),
    "); falling back to column order.",
    call. = FALSE
  )

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

#' Identify threshold or intercept rows among *reported* coefficient names.
#'
#' Reported PO thresholds are `Y<=k` or `Y<C`; reported multinomial intercepts
#' keep VGAM's `(Intercept):k`. Slope names come from design-matrix columns
#' and never take either form. An optional `stage1:`/`stage2:` prefix from
#' `coef.tsco()` is tolerated.
#'
#' @param x Character vector of reported coefficient names.
#' @return Logical vector.
#' @noRd
.tsco_is_threshold_name <- function(x) {
  x <- sub("^stage[12]:", "", x)
  grepl("^Y<", x) | grepl("^\\(Intercept\\)", x)
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

.tsco_stage_coef_terms <- function(
    fit,
    formula,
    data,
    engine,
    stage,
    kind,
    response_levels,
    reported_names = NULL) {

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

  mm <- stats::model.matrix(
    stats::delete.response(tt),
    data = data
  )

  assign <- attr(mm, "assign")
  design_idx <- which(assign > 0L)
  design_names <- colnames(mm)[design_idx]
  design_terms <- term_labels[assign[design_idx]]

  if (is.null(reported_names)) {
    reported_names <- .tsco_stage_reported_names(
      fit,
      kind = kind,
      engine = engine,
      response_levels = response_levels,
      formula = formula,
      data = data
    )
  }

  out <- rep(NA_character_, length(cf))
  names(out) <- paste0(
    stage, ":",
    .tsco_apply_reported_names(coef_names, reported_names)
  )

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

  # Use raw VGAM names for design-matrix matching. Reported category labels
  # are only for the names of the returned mapping.
  base_names <- .tsco_base_coef_name(coef_names)
  slope_idx <- which(!.tsco_is_intercept_name(base_names))
  matched <- match(base_names[slope_idx], design_names)

  if (anyNA(matched)) {
    stop(
      "Cannot align vglm coefficients with formula terms: ",
      paste(
        unique(base_names[slope_idx][is.na(matched)]),
        collapse = ", "
      ),
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
      fit = object$fit_stage1,
      formula = object$stage1_formula,
      data = object$data_stage1,
      engine = object$stage1_engine,
      stage = "stage1",
      kind = object$stage1,
      response_levels = .tsco_stage_levels(object, 1L),
      reported_names = .tsco_stage_reported_names(
        object$fit_stage1,
        kind = object$stage1,
        engine = object$stage1_engine,
        response_levels = .tsco_stage_levels(object, 1L),
        formula = object$stage1_formula,
        data = object$data_stage1,
        collapsed_label = object$lower_collapsed_label,
        cutoff_level = object$cutoff_level
      )
    ),
    .tsco_stage_coef_terms(
      fit = object$fit_stage2,
      formula = object$stage2_formula,
      data = object$data_stage2,
      engine = object$stage2_engine,
      stage = "stage2",
      kind = object$stage2,
      response_levels = .tsco_stage_levels(object, 2L),
      reported_names = .tsco_stage_reported_names(
        object$fit_stage2,
        kind = object$stage2,
        engine = object$stage2_engine,
        response_levels = .tsco_stage_levels(object, 2L),
        formula = object$stage2_formula,
        data = object$data_stage2,
        collapsed_label = object$lower_collapsed_label,
        cutoff_level = object$cutoff_level
      )
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

#' Does the installed version of `rms::orm()` accept case weights?
#'
#' `rms::orm()` gained a `weights` argument only in later releases. Feature
#' detection keeps `tsco()` working across `rms` versions instead of failing
#' inside the fitter with an "unused argument" error.
#' @noRd
.tsco_orm_supports_weights <- function() {
  if (!requireNamespace("rms", quietly = TRUE)) {
    return(FALSE)
  }

  fmls <- tryCatch(names(formals(rms::orm)), error = function(e) character())

  "weights" %in% fmls
}

#' Convergence and separation diagnostics for one fitted stage.
#'
#' Complete or quasi-complete separation drives maximum-likelihood estimates
#' and their standard errors toward infinity. The fitters usually still
#' return an object, so the resulting Wald statistics look ordinary while
#' being meaningless. This helper flags the two symptoms that matter in
#' practice -- a fit the backend reports as failed, and coefficients whose
#' magnitude and standard error are both implausibly large -- so that
#' `summary.tsco()` can say so rather than printing the numbers unannotated.
#'
#' @param fit A fitted stage model.
#' @param engine The stage fitting backend, either "orm" or "vglm".
#' @param coef_table A stage coefficient table from `.tsco_coef_table()`.
#' @param max_abs_coef,max_se Thresholds above which a coefficient is treated
#'   as separated.
#' @return A list with `converged`, `separated`, and `messages`.
#' @noRd
.tsco_stage_diagnostics <- function(
    fit,
    engine,
    coef_table = NULL,
    max_abs_coef = 10,
    max_se = 50) {

  messages <- character()
  converged <- TRUE

  if (identical(engine, "orm")) {
    if (isTRUE(fit$fail)) {
      converged <- FALSE
      messages <- c(messages, "the rms::orm() fit did not converge")
    }
  } else if (identical(engine, "vglm")) {
    iter <- tryCatch(methods::slot(fit, "iter"), error = function(e) NULL)
    maxit <- tryCatch(fit@control$maxit, error = function(e) NULL)

    if (is.numeric(iter) && is.numeric(maxit) &&
        length(iter) == 1L && length(maxit) == 1L &&
        iter >= maxit) {
      converged <- FALSE
      messages <- c(
        messages,
        paste0("the VGAM::vglm() fit reached its iteration limit (", maxit, ")")
      )
    }
  }

  separated <- FALSE

  if (!is.null(coef_table) &&
      is.matrix(coef_table) &&
      all(c("Estimate", "Std. Error") %in% colnames(coef_table))) {

    # The heuristic is applied to slope rows only. Threshold rows are
    # excluded because `rms::vcov.orm()` may legitimately omit their
    # covariance entries, which `.tsco_stage_coef_table()` pads with `NA`; a
    # nonfinite threshold standard error is therefore not evidence of
    # separation, whereas a nonfinite slope standard error is.
    is_slope <- !.tsco_is_threshold_name(rownames(coef_table))

    b <- coef_table[, "Estimate"]
    se <- coef_table[, "Std. Error"]

    flagged <- which(
      is_slope & (
        !is.finite(b) |
          !is.finite(se) |
          (abs(b) > max_abs_coef & se > max_se)
      )
    )

    if (length(flagged) > 0L) {
      separated <- TRUE
      messages <- c(
        messages,
        paste0(
          "coefficient(s) ",
          paste(rownames(coef_table)[flagged], collapse = ", "),
          " have extreme or nonfinite estimates and standard errors, a ",
          "heuristic signature of complete or quasi-complete separation"
        )
      )
    }
  }

  list(
    converged = converged,
    separated = separated,
    messages = messages
  )
}

#' Rebuild a formula right-hand side from model-frame column names.
#'
#' `tsco()` fits its two stages against data frames whose columns are the
#' *evaluated* model-frame columns, so a column produced by `log(x)` is named
#' `log(x)` and the variable `x` no longer exists. Reusing the original
#' right-hand side would therefore fail to resolve. Rebuilding each term from
#' the backticked model-frame names lets inline transformations such as
#' `log(x)`, `poly(x, 2)` and `rms::rcs(x)` work, while leaving formulas made
#' only of plain variables completely unchanged.
#'
#' @param tt A terms object for the user's original formula.
#' @return A character vector of term labels for the stage formulas.
#' @noRd
.tsco_stage_rhs_labels <- function(tt) {
  term_labels <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")

  labels <- character()

  if (length(term_labels) > 0L) {
    if (is.null(fac) || length(fac) == 0L || is.null(rownames(fac))) {
      labels <- paste0("`", term_labels, "`")
    } else {
      labels <- vapply(
        term_labels,
        function(lab) {
          vars <- rownames(fac)[fac[, lab] > 0L]
          paste0("`", vars, "`", collapse = ":")
        },
        character(1L),
        USE.NAMES = FALSE
      )
    }
  }

  # Offsets are rejected in `tsco()` before this point, so there is nothing to
  # carry through. Rebuilding them here would be the only untested route by
  # which one could reach the fitters.
  labels
}

.tsco_stage_engine <- function(kind, grouped, po_engine, weighted = FALSE) {
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

  if (identical(po_engine, "vglm")) {
    return("vglm")
  }

  # `weights` requires a backend that accepts them.
  if (isTRUE(weighted) && !.tsco_orm_supports_weights()) {
    if (identical(po_engine, "orm")) {
      stop(
        "`po_engine = \"orm\"` cannot be used with `weights` because the ",
        "installed version of rms (", utils::packageVersion("rms"), ") does ",
        "not accept case weights in `orm()`; use `po_engine = \"vglm\"` or ",
        "upgrade rms.",
        call. = FALSE
      )
    }

    return("vglm")
  }

  "orm"
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

    # `rms::orm()` warns on every weighted fit that weights are ignored by
    # `rms::validate()` and `rms::bootcov()`. `tsco()` calls neither, and the
    # caveat is stated in `?tsco` under `weights`, so the note is muffled here
    # rather than repeated twice per fit and twice per LRT refit. Only that
    # exact message is silenced; any other warning from the fitter passes.
    fit <- withCallingHandlers(
      do.call(rms::orm, c(args, stage_args)),
      warning = function(w) {
        if (grepl("weights are ignored in model validation",
                  conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )

    return(.tsco_trim_fit_call(fit, fun = quote(rms::orm)))
  }

  if (!identical(engine, "vglm")) {
    stop("Unknown stage fitting engine: ", engine, call. = FALSE)
  }

  if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop("Package 'VGAM' is required for this stage.", call. = FALSE)
  }

  family_call <- if (identical(kind, "po")) {
    quote(VGAM::cumulative(link = "logitlink", parallel = TRUE, reverse = FALSE))
  } else if (identical(kind, "multinomial")) {
    quote(VGAM::multinomial(refLevel = 1))
  } else {
    stop("Unknown stage model kind: ", kind, call. = FALSE)
  }

  family <- eval(family_call)

  args <- list(
    formula = formula,
    family = family,
    data = data
  )

  if (!is.null(weights)) {
    args$weights <- weights
  }

  .tsco_trim_fit_call(
    do.call(VGAM::vglm, c(args, stage_args)),
    fun = quote(VGAM::vglm),
    family = family_call
  )
}

#' Replace bulky literals in a fitted stage's recorded call.
#'
#' `do.call()` splices the stage data frame, the case weights and the family
#' object into the call each fitter records, which inflates every fitted
#' `tsco` object by roughly the size of the analysis data. Nothing in the
#' package re-evaluates those recorded calls -- likelihood-ratio refits build
#' fresh calls and prediction goes through `newdata` -- so the literals are
#' replaced with placeholder symbols after fitting.
#'
#' @param fit A fitted stage model.
#' @param fun A call or name to record as the fitted call's function.
#' @param family A replacement expression for the recorded `family` argument.
#' @noRd
.tsco_trim_fit_call <- function(fit, fun = NULL, family = NULL) {
  trim <- function(cl) {
    if (!is.call(cl)) {
      return(cl)
    }

    # `do.call()` also splices in the fitter function itself, so the recorded
    # call carries a copy of its whole body.
    if (!is.null(fun)) {
      cl[[1L]] <- fun
    }

    for (nm in intersect(c("data", "weights"), names(cl))) {
      cl[[nm]] <- as.name(paste0(".tsco_", nm))
    }

    if (!is.null(family) && "family" %in% names(cl)) {
      cl[["family"]] <- family
    }

    cl
  }

  if (isS4(fit)) {
    if (methods::.hasSlot(fit, "call")) {
      methods::slot(fit, "call") <- trim(methods::slot(fit, "call"))
    }

    return(fit)
  }

  if (!is.null(fit$call)) {
    fit$call <- trim(fit$call)
  }

  fit
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

  # `stats::drop.terms()` cannot build an intercept-only formula on every
  # supported R version, so construct the reduced formula directly. This keeps
  # joint tests available for single-predictor models such as `y ~ treatment`.
  if (length(term_labels) == 1L) {
    if (attr(tt, "intercept") == 0L) {
      stop(
        "Cannot drop the only term from a formula without an intercept.",
        call. = FALSE
      )
    }

    f_reduced <- stats::as.formula(
      paste0(.tsco_deparse(formula[[2L]]), " ~ 1"),
      env = environment(formula)
    )

    return(f_reduced)
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

#' Log-likelihood and degrees of freedom of an intercept-only stage.
#'
#' A PO or multinomial model with no covariates is saturated in the category
#' probabilities, so its maximum-likelihood fit is the observed (weighted)
#' category distribution and its log-likelihood has a closed form. Using that
#' form avoids refitting an intercept-only model, which not every backend
#' supports: `rms::orm(y ~ 1)` fails outright on some versions of rms. The
#' value omits multinomial combinatorial constants, matching
#' `.tsco_count_loglik()` and therefore `object$logLik_stage1` and
#' `object$logLik_stage2`. It is *not* comparable with a grouped-count
#' `VGAM::vglm()` fit's own `@criterion["loglikelihood"]`, which includes
#' those constants; `.tsco_joint_lrt_tests()` never uses backend
#' log-likelihoods for that reason.
#'
#' @param totals Weighted category totals for the stage response.
#' @return A list with `loglik` and `df`.
#' @noRd
.tsco_null_stage_fit_stats <- function(totals) {
  totals <- as.numeric(totals)

  if (length(totals) < 2L || any(!is.finite(totals)) || sum(totals) <= 0) {
    return(list(loglik = NA_real_, df = NA_real_))
  }

  p <- totals / sum(totals)
  positive <- totals > 0

  list(
    loglik = sum(totals[positive] * log(p[positive])),
    df = length(totals) - 1L
  )
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

  # Every log-likelihood entering a likelihood-ratio statistic is computed on
  # the package's constant-free scale (`.tsco_count_loglik()`), never taken
  # from a backend. The backends do not agree with each other or with the
  # closed-form intercept-only fit: `VGAM::vglm()` on a grouped-count response
  # includes the multinomial combinatorial constants in `@criterion`, while
  # `rms::orm()` and the closed form omit them. Those constants cancel between
  # two refits of the same grouped data, but not between a grouped `vglm` fit
  # and the closed form, so mixing scales made a grouped `y ~ treatment` LRT
  # wildly anti-conservative.
  ll1_full <- object$logLik_stage1
  ll2_full <- object$logLik_stage2

  df1_full <- .tsco_df(object$fit_stage1)
  df2_full <- .tsco_df(object$fit_stage2)

  if (is.null(ll1_full) || is.null(ll2_full) ||
      !is.finite(ll1_full) || !is.finite(ll2_full)) {
    stop(
      "Could not extract full-model log-likelihoods for LRTs.",
      call. = FALSE
    )
  }

  if (is.null(object$counts_stage1) || is.null(object$counts_stage2)) {
    stop(
      "Joint LRTs require the stage count matrices stored in the tsco object.",
      call. = FALSE
    )
  }

  stage1_levels <- .tsco_stage_levels(object, 1L)
  stage2_levels <- .tsco_stage_levels(object, 2L)

  # Score a reduced stage fit against the stored weighted counts.
  reduced_loglik <- function(fit, engine, data, counts, levels) {
    p <- .tsco_predict_stage_prob(
      fit,
      engine = engine,
      newdata = data,
      levels = levels
    )
    p <- .tsco_align_prob(p, levels)
    .tsco_count_loglik(counts, p)
  }

  stage1_args <- object$stage1_args
  stage2_args <- object$stage2_args

  # Dropping the only term leaves an intercept-only model in both stages,
  # whose fit is available in closed form.
  null_reduced <- length(term_labels) == 1L

  res <- lapply(term_labels, function(term) {
    tryCatch({
      if (null_reduced) {
        null1 <- .tsco_null_stage_fit_stats(object$totals_stage1)
        null2 <- .tsco_null_stage_fit_stats(object$totals_stage2)

        ll1_reduced <- null1$loglik
        ll2_reduced <- null2$loglik
        df1_reduced <- null1$df
        df2_reduced <- null2$df

        if (!is.finite(ll1_reduced) || !is.finite(ll2_reduced)) {
          stop(
            "Could not compute intercept-only reduced-model log-likelihoods.",
            call. = FALSE
          )
        }
      } else {

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

      ll1_reduced <- reduced_loglik(
        fit1_reduced, stage1_engine, object$data_stage1,
        object$counts_stage1, stage1_levels
      )
      ll2_reduced <- reduced_loglik(
        fit2_reduced, stage2_engine, object$data_stage2,
        object$counts_stage2, stage2_levels
      )

      if (!is.finite(ll1_reduced) || !is.finite(ll2_reduced)) {
        stop("Could not compute reduced-model log-likelihoods.", call. = FALSE)
      }

      df1_reduced <- .tsco_df(fit1_reduced)
      df2_reduced <- .tsco_df(fit2_reduced)

      }

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

  # If `newdata` supplies the raw variables rather than the evaluated
  # model-frame columns, evaluate the original right-hand side against it so
  # that inline transformations such as `log(x)` resolve the same way they did
  # at fitting time. When `newdata` already holds the model-frame columns
  # (including `object$predict_data` itself) it is used as supplied.
  mf_cols <- names(object$predict_data)

  if (!is.null(object$terms_rhs) &&
      length(mf_cols) > 0L &&
      !all(mf_cols %in% names(nd))) {

    nd <- tryCatch(
      stats::model.frame(
        object$terms_rhs,
        data = nd,
        na.action = stats::na.pass
      ),
      error = function(e) {
        stop(
          "Could not evaluate the model formula on `newdata`: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

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

#' Rows of prediction data with a factor level absent from the fitted data.
#'
#' `.tsco_prepare_newdata()` rejects such levels in user-supplied `newdata`,
#' so this only bites for `predict(fit)` without `newdata`, where
#' `predict_data` still holds rows that were removed before fitting because
#' their weight was exactly zero. Neither stage carries information about a
#' level that occurs only in those rows, so nothing about them is
#' identifiable and `predict()` returns `NA` for the whole row.
#'
#' @param object A fitted `tsco` object.
#' @param newdata Prepared prediction data.
#' @return A logical vector, `TRUE` for rows that cannot be predicted.
#' @noRd
.tsco_unfitted_rows <- function(object, newdata) {
  observed <- object$predictor_observed_levels
  n <- nrow(newdata)

  if (is.null(observed) || length(observed) == 0L) {
    return(rep(FALSE, n))
  }

  out <- rep(FALSE, n)

  for (nm in names(observed)) {
    if (!nm %in% names(newdata)) {
      next
    }

    vals <- as.character(newdata[[nm]])
    out <- out | (!is.na(vals) & !(vals %in% observed[[nm]]))
  }

  out
}

#' Rows of prediction data with a missing value in any model predictor.
#'
#' Only the model-frame predictor columns count; extra columns in `newdata`
#' are ignored. An intercept-only model has no predictors and no incomplete
#' rows.
#'
#' @param object A fitted `tsco` object.
#' @param newdata Prepared prediction data.
#' @return A logical vector, `TRUE` for rows that cannot be predicted.
#' @noRd
.tsco_incomplete_rows <- function(object, newdata) {
  cols <- intersect(names(object$predict_data), names(newdata))

  if (length(cols) == 0L) {
    return(rep(FALSE, nrow(newdata)))
  }

  .tsco_bad_cells(newdata[, cols, drop = FALSE])$rows
}

#' Locate missing or non-finite predictor values in a model frame.
#'
#' A numeric or matrix column (`poly()`, `splines::ns()`) is bad where it is
#' `NA`, `NaN` or infinite; any other column (factor, character, logical) is
#' bad where it is `NA`. Infinite values arise from transformations such as
#' `log(0)` and are just as unusable as `NA`: a fitter fails on them, and a
#' prediction at an infinite linear predictor is a degenerate 0/1 row, not a
#' model-based probability.
#'
#' @param df A data frame of predictor columns (a model frame minus the
#'   response, possibly with matrix-valued columns).
#' @return A list with `rows`, a logical vector marking rows with any bad
#'   cell, and `columns`, the names of the columns containing one.
#' @noRd
.tsco_bad_cells <- function(df) {
  n <- nrow(df)
  rows <- rep(FALSE, n)
  columns <- character()

  for (nm in names(df)) {
    x <- df[[nm]]

    bad <- if (is.numeric(x) || is.complex(x)) {
      !is.finite(x)
    } else {
      is.na(x)
    }

    if (is.matrix(bad)) {
      bad <- rowSums(bad) > 0
    }

    bad <- as.logical(bad)
    bad[is.na(bad)] <- TRUE

    if (any(bad)) {
      columns <- c(columns, nm)
      rows <- rows | bad
    }
  }

  list(rows = rows, columns = columns)
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


#' Empirical distribution of the lower categories among stage-1 observations.
#'
#' Used by `predict(unsupported_stage1 = "empirical")`. The distribution is
#' taken from the weighted category totals stored at fitting time, so case
#' weights are honoured and, consistent with the rest of the package, a
#' declared lower level with no observations receives probability exactly
#' zero rather than a share of a uniform split.
#'
#' @param object A fitted `tsco` object.
#' @return A named probability vector over `object$lower_levels`.
#' @noRd
.tsco_empirical_lower_probs <- function(object) {
  totals <- object$totals_stage1

  if (is.null(totals) || is.null(names(totals))) {
    stop(
      "Cannot compute empirical lower-category probabilities because ",
      "stage-1 category totals are not stored in the TSCO object.",
      call. = FALSE
    )
  }

  lower <- object$lower_levels

  counts <- stats::setNames(numeric(length(lower)), lower)
  common <- intersect(names(totals), lower)
  counts[common] <- as.numeric(totals[common])

  if (any(!is.finite(counts)) || sum(counts) <= 0) {
    warning(
      "Could not compute empirical lower-category probabilities; ",
      "using a uniform split among the observed lower categories instead.",
      call. = FALSE
    )

    observed <- intersect(.tsco_stage_levels(object, 1L), lower)
    counts[] <- 0
    counts[observed] <- 1
  }

  counts / sum(counts)
}


#' Design-matrix column names for a stage, excluding the intercept.
#'
#' Used to give proportional-odds slopes the same reported names regardless of
#' whether the stage was fitted with `rms::orm()` or `VGAM::vglm()`.
#'
#' @param formula The stage formula.
#' @param data The stage fitting data.
#' @noRd
.tsco_stage_design_names <- function(formula, data) {
  if (is.null(formula) || is.null(data)) {
    return(NULL)
  }

  tryCatch({
    tt <- stats::delete.response(stats::terms(formula))
    mm <- stats::model.matrix(tt, data = data)
    assign <- attr(mm, "assign")
    colnames(mm)[assign > 0L]
  }, error = function(e) NULL)
}

#' Reported coefficient names for a proportional-odds stage.
#'
#' `coef.tsco()` reports PO thresholds on the lower-tail scale
#' `Pr(Y <= level_j | X) = expit(gamma_j + x'beta)`, so the thresholds are
#' labelled by the level they bound rather than by the backend's own
#' convention. `rms::orm()` names its thresholds for the *upper* tail
#' (`y>=2`), which contradicts the reported value once the sign correction in
#' `.tsco_sign_vector()` has been applied. Slopes take the design-matrix
#' column names, which `VGAM::vglm()` already uses, so that both engines
#' produce identical names for the same model.
#'
#' @param coef_names The backend's coefficient names.
#' @param response_levels The levels of the stage's response variable.
#' @param n_thresholds The number of threshold (intercept) coefficients, which
#'   come first in both backends' coefficient vectors.
#' @param design_names Design-matrix column names for the slopes, excluding
#'   the intercept. Optional; raw names are kept when absent.
#' @param collapsed_label,cutoff_level The stage-2 collapsed level's internal
#'   label and the cutoff level it stands for, so that its threshold can be
#'   reported as `Y<C`. Optional.
#' @noRd
.tsco_po_reported_coef_names <- function(
    coef_names,
    response_levels,
    n_thresholds,
    design_names = NULL,
    collapsed_label = NULL,
    cutoff_level = NULL) {

  out <- coef_names

  # In stage 2 the first level is the collapsed lower partition, so the
  # threshold that bounds it is exactly `Y < C`; naming it `Y<=<internal
  # label>` would be both ugly and confusing.
  threshold_labels <- paste0("Y<=", response_levels)

  if (!is.null(collapsed_label) && !is.null(cutoff_level)) {
    is_collapsed <- response_levels == collapsed_label
    threshold_labels[is_collapsed] <- paste0("Y<", cutoff_level)
  }

  if (!is.finite(n_thresholds) || n_thresholds < 1L) {
    return(out)
  }

  n_thresholds <- as.integer(n_thresholds)

  if (n_thresholds > length(coef_names) ||
      n_thresholds >= length(response_levels)) {
    return(out)
  }

  out[seq_len(n_thresholds)] <- threshold_labels[seq_len(n_thresholds)]

  slope_idx <- seq.int(n_thresholds + 1L, length(coef_names))

  if (!is.null(design_names) && length(design_names) == length(slope_idx)) {
    out[slope_idx] <- design_names
  }

  out
}

#' Reported coefficient names for one fitted TSCO stage.
#'
#' Returns a character vector of reported names, named by the backend's own
#' coefficient names so that callers can look up names for an arbitrary
#' subset or ordering of coefficients.
#'
#' @param fit A fitted stage model.
#' @param kind "po" or "multinomial".
#' @param engine The stage fitting backend, either "orm" or "vglm".
#' @param response_levels The levels of the stage's response variable.
#' @param formula,data The stage formula and fitting data, used to recover
#'   design-matrix column names for PO slopes. Optional.
#' @param collapsed_label,cutoff_level The stage-2 collapsed level's internal
#'   label and the cutoff level it stands for. Optional.
#' @noRd
.tsco_stage_reported_names <- function(
    fit,
    kind,
    engine = NULL,
    response_levels = NULL,
    formula = NULL,
    data = NULL,
    collapsed_label = NULL,
    cutoff_level = NULL) {

  coef_names <- names(stats::coef(fit))

  if (is.null(coef_names)) {
    return(NULL)
  }

  if (is.null(engine)) {
    engine <- if (inherits(fit, "orm")) "orm" else "vglm"
  }

  if (is.null(response_levels)) {
    out <- coef_names
    names(out) <- coef_names
    return(out)
  }

  if (identical(kind, "multinomial")) {
    out <- .tsco_reported_coef_names(
      coef_names,
      kind = kind,
      response_levels = response_levels
    )

    names(out) <- coef_names
    return(out)
  }

  n_thresholds <- if (identical(engine, "orm")) {
    fit$non.slopes
  } else {
    sum(grepl("(Intercept)", coef_names, fixed = TRUE))
  }

  out <- .tsco_po_reported_coef_names(
    coef_names,
    response_levels = response_levels,
    n_thresholds = n_thresholds,
    design_names = .tsco_stage_design_names(formula, data),
    collapsed_label = collapsed_label,
    cutoff_level = cutoff_level
  )

  names(out) <- coef_names
  out
}

.tsco_reported_coef_names <- function(
    coef_names,
    kind,
    response_levels) {

  if (!identical(kind, "multinomial")) {
    return(coef_names)
  }

  if (length(response_levels) < 2L) {
    stop(
      "Internal error: multinomial response must have at least two levels.",
      call. = FALSE
    )
  }

  # refLevel = 1, so linear predictors correspond to response_levels[-1].
  target_levels <- response_levels[-1L]

  has_suffix <- grepl(":[0-9]+$", coef_names)

  # VGAM may omit the :1 suffix for binary multinomial models.
  if (!any(has_suffix)) {
    if (length(target_levels) == 1L) {
      return(paste0(coef_names, ":", target_levels))
    }

    stop(
      "Internal error: multinomial coefficient names lack linear-predictor ",
      "suffixes for a response with multiple nonreference categories.",
      call. = FALSE
    )
  }

  index <- rep(NA_integer_, length(coef_names))

  index[has_suffix] <- as.integer(
    sub("^.*:([0-9]+)$", "\\1", coef_names[has_suffix])
  )

  bad <- has_suffix & (
    is.na(index) |
      index < 1L |
      index > length(target_levels)
  )

  if (any(bad)) {
    stop(
      "Internal error: multinomial coefficient suffixes could not be ",
      "matched to response levels.",
      call. = FALSE
    )
  }

  # In multicategory VGAM fits, all coefficient names should ordinarily carry
  # a predictor suffix. Treat mixed suffixed/unsuffixed names as an internal
  # inconsistency rather than silently returning ambiguous names.
  if (any(!has_suffix)) {
    stop(
      "Internal error: multinomial coefficient names contain a mixture of ",
      "suffixed and unsuffixed values.",
      call. = FALSE
    )
  }

  paste0(
    sub(":[0-9]+$", "", coef_names),
    ":",
    target_levels[index]
  )
}

# -------------------------------------------------------------------------
# Declared-but-unobserved outcome levels
# -------------------------------------------------------------------------

#' Which of a stage's declared levels actually carry observations.
#'
#' @param y A stage response: a factor, or a matrix of category counts.
#' @param levels The stage's declared levels, in order.
#' @noRd
.tsco_observed_levels <- function(y, levels) {
  if (is.matrix(y) || is.data.frame(y)) {
    totals <- colSums(as.matrix(y))
    return(levels[levels %in% colnames(as.matrix(y))[totals > 0]])
  }

  counts <- tabulate(match(as.character(y), levels), nbins = length(levels))
  levels[counts > 0]
}

#' Insert zero-probability columns for declared levels that were not fitted.
#'
#' A level declared in `levels` but absent from the data has a maximum
#' likelihood fitted probability of exactly zero, and its threshold sits on the
#' boundary of the parameter space. Neither backend keeps such a level, so the
#' fitted probability matrix spans only the observed levels; this restores the
#' declared layout so that predictions always have one column per outcome level
#' and rows that still sum to one.
#'
#' @param p Probability matrix over `observed`, in that column order.
#' @param observed,declared Character vectors of level labels.
#' @noRd
.tsco_expand_prob <- function(p, observed, declared) {
  p <- as.matrix(p)

  if (identical(observed, declared)) {
    return(p)
  }

  out <- matrix(
    0,
    nrow = nrow(p),
    ncol = length(declared),
    dimnames = list(rownames(p), declared)
  )

  out[, observed] <- p[, observed, drop = FALSE]
  out
}

#' Drop declared-but-unobserved levels from a stage response before fitting.
#'
#' `rms::orm()` fails with an opaque `'names' attribute` error, and
#' `VGAM::vglm()` silently drops the level and returns a probability matrix of
#' the wrong width. Neither is recoverable downstream, so the level is removed
#' here, recorded, and restored at prediction time by `.tsco_expand_prob()`.
#'
#' @param y A stage response: a factor, or a matrix of category counts.
#' @param levels The stage's declared levels, in order.
#' @param stage_label Text naming the stage, used in messages.
#' @return A list with the reduced `y`, the `observed` levels, and the
#'   `dropped` ones.
#' @noRd
.tsco_drop_unobserved_levels <- function(y, levels, stage_label) {
  observed <- .tsco_observed_levels(y, levels)
  dropped <- setdiff(levels, observed)

  if (length(observed) < 2L) {
    stop(
      stage_label, " has ", length(observed),
      " outcome level(s) with observations (",
      if (length(observed) == 0L) "none" else paste(observed, collapse = ", "),
      ") out of the ", length(levels), " declared (",
      paste(levels, collapse = ", "),
      "). At least two are needed to fit a model. Supply more data or a ",
      "different `cutoff_level`.",
      call. = FALSE
    )
  }

  if (length(dropped) == 0L) {
    return(list(y = y, observed = observed, dropped = character()))
  }

  warning(
    stage_label, ": outcome level(s) ", paste(dropped, collapse = ", "),
    " are declared in `levels` but have no observations. They are dropped ",
    "from the fit and carried through predictions with probability zero; ",
    "their thresholds are at the boundary of the parameter space. Estimates ",
    "and tests are conditional on the observed outcome levels.",
    call. = FALSE
  )

  if (is.matrix(y) || is.data.frame(y)) {
    y <- as.matrix(y)[, observed, drop = FALSE]
  } else {
    y <- factor(
      as.character(y),
      levels = observed,
      ordered = is.ordered(y)
    )
  }

  list(y = y, observed = observed, dropped = dropped)
}

#' A stage's observed levels, falling back to its declared levels.
#'
#' @param object A "tsco" object.
#' @param stage 1 or 2.
#' @noRd
.tsco_stage_levels <- function(object, stage) {
  declared <- if (stage == 1L) object$lower_levels else object$collapsed_levels
  observed <- if (stage == 1L) {
    object$stage1_observed_levels
  } else {
    object$stage2_observed_levels
  }

  if (is.null(observed)) declared else observed
}
