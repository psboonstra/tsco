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

.tsco_coef_table <- function(fit) {
  s <- tryCatch(summary(fit), error = function(e) NULL)

  if (!is.null(s)) {
    ct <- tryCatch(stats::coef(s), error = function(e) NULL)

    if (is.matrix(ct) || is.data.frame(ct)) {
      return(as.matrix(ct))
    }

    # VGAM summary objects are commonly S4. Try common slots if coef() failed.
    if (methods::is(s, "summaryvglm")) {
      sn <- methods::slotNames(s)

      for (slot_name in c("coef3", "coef4", "coeftable")) {
        if (slot_name %in% sn) {
          ct <- tryCatch(methods::slot(s, slot_name), error = function(e) NULL)

          if (is.matrix(ct) || is.data.frame(ct)) {
            return(as.matrix(ct))
          }
        }
      }
    }
  }

  # Fallback: estimates only.
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)

  if (is.null(cf)) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  matrix(
    cf,
    ncol = 1L,
    dimnames = list(names(cf), "Estimate")
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
