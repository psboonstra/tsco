# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.tsco_deparse <- function(x, width.cutoff = 500L) {
  paste(deparse(x, width.cutoff = width.cutoff), collapse = " ")
}

.tsco_df <- function(fit) {
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)

  if (is.null(cf)) {
    return(NA_real_)
  }

  length(cf)
}


.tsco_stage_coef_table <- function(fit, kind, po.reverse = FALSE) {
  b_raw <- tryCatch(stats::coef(fit), error = function(e) NULL)
  V_raw <- tryCatch(stats::vcov(fit), error = function(e) NULL)

  if (is.null(b_raw)) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  s <- .tsco_sign_vector(fit, kind = kind, po.reverse = po.reverse)

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

#' Construct a stage-level coefficient table on the TSCO/manuscript scale.
#'
#' Estimates and covariance matrices are transformed using `.tsco_sign_vector()`.
#' Standard errors and p-values are computed from the transformed covariance
#' matrix.
#'
#' @param fit A fitted stage model (a "vglm" object).
#' @param kind "po" or "multinomial", i.e. object$stage1 or object$stage2.
#' @param po.reverse Logical. Whether PO stages were fit with
#'   `VGAM::cumulative(reverse = TRUE)`. If `TRUE`, non-intercept PO
#'   coefficients are sign-flipped to match the manuscript's beta convention.
.tsco_coef_table <- function(fit, kind, po.reverse = FALSE) {
  .tsco_stage_coef_table(
    fit = fit,
    kind = kind,
    po.reverse = po.reverse
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

#' Sign vector correcting a fitted stage's coefficients to the manuscript's
#' parameterization.
#'
#' The manuscript parameterizes both submodels with a subtractive
#' convention on x'beta (Eq. 1 for PO, Eq. 2 for MR/multinomial). As
#' actually called by tsco()'s make_family():
#'   - "po" (VGAM::cumulative(reverse = po.reverse)): with
#'     po.reverse = FALSE (the package default as of this fix), VGAM's
#'     Pr(Y <= j) = expit(alpha_j + x'beta) is algebraically equivalent to
#'     Pr(Y >= j+1) = expit(-alpha_j - x'beta), which already matches the
#'     manuscript's sign on x'beta once cutpoints are relabeled. No flip
#'     needed. (Confirmed by simulation in test-sign_recovery_simulation.R;
#'     if po.reverse is ever changed back to TRUE, this would need to
#'     change too.)
#'   - "multinomial" (VGAM::multinomial(refLevel = 1)): uses the standard
#'     additive multinomial-logit convention, eta = alpha + x'beta, which
#'     is the *opposite* sign of the manuscript's Eq. 2. Confirmed by
#'     simulation in test-multinomial_sign_recovery.R (recovered
#'     coefficient was the negative of the true simulated value). Every
#'     non-intercept coefficient needs flipping.
#'
#' This helper is the single point of control for that correction so that,
#' if make_family()'s arguments ever change, there is one place to update
#' the sign logic rather than several call sites (coef.tsco, vcov.tsco,
#' summary.tsco).
#'
#' @param fit A fitted stage model (a "vglm" object).
#' @param kind "po" or "multinomial", i.e. object$stage1 or object$stage2.
#' @param po.reverse Logical. Whether PO stages were fit with
#'   `VGAM::cumulative(reverse = TRUE)`. If `TRUE`, non-intercept PO
#'   coefficients are sign-flipped to match the manuscript's beta convention.
#' @return A named numeric vector of +1/-1, same length and names as
#'   stats::coef(fit).
.tsco_sign_vector <- function(fit, kind, po.reverse = FALSE) {
  cf <- stats::coef(fit)
  s <- rep(1, length(cf))
  names(s) <- names(cf)

  is_intercept <- grepl("(Intercept)", names(cf), fixed = TRUE)

  if (identical(kind, "multinomial")) {
    # VGAM multinomial uses alpha + x'beta, while manuscript uses
    # alpha - x'beta. Flip non-intercepts.
    s[!is_intercept] <- -1
  }

  if (identical(kind, "po") && isTRUE(po.reverse)) {
    # With reverse = TRUE, VGAM's PO slopes have the opposite sign from
    # the manuscript beta convention. Flip non-intercepts.
    s[!is_intercept] <- -1
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
  p <- as.matrix(p)

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

  base_names <- .tsco_base_coef_name(names(b))

  test_names <- unique(base_names[!.tsco_is_intercept_name(base_names)])

  if (length(test_names) == 0L) {
    return(data.frame())
  }

  res <- lapply(test_names, function(term) {
    idx <- which(base_names == term)

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

.tsco_make_family <- function(kind, po.reverse = FALSE) {
  if (identical(kind, "po")) {
    VGAM::cumulative(
      link = "logitlink",
      parallel = TRUE,
      reverse = po.reverse
    )
  } else if (identical(kind, "multinomial")) {
    VGAM::multinomial(refLevel = 1)
  } else {
    stop("Unknown stage model kind: ", kind, call. = FALSE)
  }
}

.tsco_fit_vglm_internal <- function(
    formula,
    family,
    data,
    weights = NULL,
    vglm_args = list()) {

  args <- list(
    formula = formula,
    family = family,
    data = data
  )

  if (!is.null(weights)) {
    args$weights <- weights
  }

  args <- c(args, vglm_args)

  do.call(VGAM::vglm, args)
}

.tsco_vglm_loglik <- function(fit) {
  out <- tryCatch({
    as.numeric(fit@criterion["loglikelihood"])
  }, error = function(e) NA_real_)

  if (length(out) != 1L || !is.finite(out)) {
    NA_real_
  } else {
    out
  }
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

.tsco_joint_lrt_tests <- function(object) {
  if (is.null(object$data_stage1) ||
      is.null(object$data_stage2) ||
      is.null(object$weights_stage1) && !("weights_stage1" %in% names(object)) ||
      is.null(object$weights_stage2) && !("weights_stage2" %in% names(object))) {
    stop(
      "Joint LRTs require stage-specific data stored in the tsco object. ",
      "Please refit the model with the current version of tsco().",
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

  ll1_full <- .tsco_vglm_loglik(object$fit_stage1)
  ll2_full <- .tsco_vglm_loglik(object$fit_stage2)

  df1_full <- .tsco_df(object$fit_stage1)
  df2_full <- .tsco_df(object$fit_stage2)

  if (!is.finite(ll1_full) || !is.finite(ll2_full)) {
    stop(
      "Could not extract full-model log-likelihoods for LRTs.",
      call. = FALSE
    )
  }

  vglm_args <- object$vglm_args
  if (is.null(vglm_args)) {
    vglm_args <- list()
  }

  fam1 <- .tsco_make_family(object$stage1, po.reverse = object$po.reverse)
  fam2 <- .tsco_make_family(object$stage2, po.reverse = object$po.reverse)

  res <- lapply(term_labels, function(term) {
    out <- tryCatch({
      f1_reduced <- .tsco_drop_term_formula(object$stage1_formula, term)
      f2_reduced <- .tsco_drop_term_formula(object$stage2_formula, term)

      fit1_reduced <- .tsco_fit_vglm_internal(
        formula = f1_reduced,
        family = fam1,
        data = object$data_stage1,
        weights = object$weights_stage1,
        vglm_args = vglm_args
      )

      fit2_reduced <- .tsco_fit_vglm_internal(
        formula = f2_reduced,
        family = fam2,
        data = object$data_stage2,
        weights = object$weights_stage2,
        vglm_args = vglm_args
      )

      ll1_reduced <- .tsco_vglm_loglik(fit1_reduced)
      ll2_reduced <- .tsco_vglm_loglik(fit2_reduced)

      df1_reduced <- .tsco_df(fit1_reduced)
      df2_reduced <- .tsco_df(fit2_reduced)

      lrt <- 2 * (
        (ll1_full + ll2_full) -
          (ll1_reduced + ll2_reduced)
      )

      # Numerical cleanup for tiny negative values from optimizer tolerance.
      if (is.finite(lrt) && lrt < 0 && lrt > -1e-7) {
        lrt <- 0
      }

      df <- (df1_full - df1_reduced) + (df2_full - df2_reduced)

      p <- if (is.finite(lrt) && is.finite(df) && df > 0) {
        stats::pchisq(lrt, df = df, lower.tail = FALSE)
      } else {
        NA_real_
      }

      data.frame(
        term = term,
        df = df,
        Chisq = lrt,
        `Pr(>Chisq)` = p,
        check.names = FALSE
      )
    }, error = function(e) {
      data.frame(
        term = term,
        df = NA_real_,
        Chisq = NA_real_,
        `Pr(>Chisq)` = NA_real_,
        check.names = FALSE
      )
    })

    out
  })

  out <- do.call(rbind, res)
  rownames(out) <- out$term
  out$term <- NULL

  out
}
