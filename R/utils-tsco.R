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

#' Extract a stage's coefficient table, sign-corrected to match
#' coef.tsco()/vcov.tsco().
#'
#' `Estimate` and any test-statistic column built by dividing the estimate
#' by its (sign-invariant) standard error -- e.g. `z value`, `t value` --
#' are flipped for non-intercept rows when `kind == "multinomial"`.
#' `Std. Error` and p-value columns are unaffected, since neither depends
#' on the sign of the estimate.
#'
#' @param fit A fitted stage model (a "vglm" object).
#' @param kind "po" or "multinomial", i.e. object$stage1 or object$stage2.
.tsco_coef_table <- function(fit, kind) {
  ct <- NULL

  s <- tryCatch(summary(fit), error = function(e) NULL)

  if (!is.null(s)) {
    ct <- tryCatch(stats::coef(s), error = function(e) NULL)

    if (!(is.matrix(ct) || is.data.frame(ct))) {
      ct <- NULL

      # VGAM summary objects are commonly S4. Try common slots if coef() failed.
      if (methods::is(s, "summaryvglm")) {
        sn <- methods::slotNames(s)

        for (slot_name in c("coef3", "coef4", "coeftable")) {
          if (slot_name %in% sn) {
            cand <- tryCatch(methods::slot(s, slot_name), error = function(e) NULL)

            if (is.matrix(cand) || is.data.frame(cand)) {
              ct <- cand
              break
            }
          }
        }
      }
    }
  }

  if (is.null(ct)) {
    # Fallback: estimates only.
    cf <- tryCatch(stats::coef(fit), error = function(e) NULL)

    if (is.null(cf)) {
      return(matrix(numeric(0), nrow = 0L, ncol = 0L))
    }

    ct <- matrix(
      cf,
      ncol = 1L,
      dimnames = list(names(cf), "Estimate")
    )
  } else {
    ct <- as.matrix(ct)
  }

  if (identical(kind, "multinomial") && nrow(ct) > 0L) {
    is_intercept <- grepl("(Intercept)", rownames(ct), fixed = TRUE)
    flip_rows <- !is_intercept

    if (any(flip_rows)) {
      cn <- colnames(ct)
      # "Estimate" and any column that is estimate-scaled (test statistics
      # such as "z value"/"t value") flip sign; "Std. Error" and p-value
      # columns do not.
      flip_cols <- grepl(
        "^Estimate$|z value|t value|Wald",
        cn,
        ignore.case = TRUE
      )

      ct[flip_rows, flip_cols] <- -ct[flip_rows, flip_cols]
    }
  }

  ct
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
#' @return A named numeric vector of +1/-1, same length and names as
#'   stats::coef(fit).
.tsco_sign_vector <- function(fit, kind) {
  cf <- stats::coef(fit)
  s <- rep(1, length(cf))
  names(s) <- names(cf)

  if (identical(kind, "multinomial")) {
    is_intercept <- grepl("(Intercept)", names(cf), fixed = TRUE)
    s[!is_intercept] <- -1
  }

  s
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
