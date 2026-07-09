# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.tsco_deparse <- function(x, width.cutoff = 500L) {
  paste(deparse(x, width.cutoff = width.cutoff), collapse = " ")
}

.tsco_safe_logLik <- function(fit) {
  out <- tryCatch({
    crit <- fit@criterion
    as.numeric(crit["loglikelihood"])
  }, error = function(e) NA_real_)

  if (length(out) == 0L || !is.finite(out)) NA_real_ else out
}

.tsco_logLik_numeric <- function(fit) {
  .tsco_safe_logLik(fit)
}

.tsco_df <- function(fit) {
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)

  if (is.null(cf)) {
    return(NA_real_)
  }

  length(cf)
}

.tsco_safe_AIC <- function(fit) {
  out <- tryCatch(VGAM::AICvlm(fit), error = function(e) NA_real_)
  if (length(out) == 0L || !is.finite(out)) NA_real_ else as.numeric(out)
}

.tsco_safe_BIC <- function(fit) {
  out <- tryCatch(VGAM::BICvlm(fit), error = function(e) NA_real_)
  if (length(out) == 0L || !is.finite(out)) NA_real_ else as.numeric(out)
}

.tsco_coef_table <- function(fit) {
  s <- tryCatch(summary(fit), error = function(e) NULL)

  if (!is.null(s)) {
    ct <- tryCatch(coef(s), error = function(e) NULL)

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
  cf <- tryCatch(coef(fit), error = function(e) NULL)

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
    stop("Internal error: response and probability matrices have different numbers of rows.")
  }

  if (ncol(y) != ncol(p)) {
    stop("Internal error: response and probability matrices have different numbers of columns.")
  }

  p <- pmax(p, eps)
  p <- p / rowSums(p)

  sum(y * log(p))
}

.tsco_factor_to_counts <- function(y, levels) {
  y <- factor(as.character(y), levels = levels)

  mm <- stats::model.matrix(~ y - 1)
  colnames(mm) <- sub("^y", "", colnames(mm))

  # model.matrix can make syntactic names, so safer to construct manually.
  out <- matrix(
    0,
    nrow = length(y),
    ncol = length(levels),
    dimnames = list(NULL, levels)
  )

  out[cbind(seq_along(y), match(as.character(y), levels))] <- 1
  out
}

.tsco_align_prob <- function(p, levels) {
  p <- as.matrix(p)

  if (is.null(colnames(p))) {
    if (ncol(p) == length(levels)) {
      colnames(p) <- levels
    } else {
      stop("Predicted probability matrix has no column names and incompatible dimension.")
    }
  }

  missing <- setdiff(levels, colnames(p))

  if (length(missing) > 0L) {
    if (ncol(p) == length(levels)) {
      colnames(p) <- levels
    } else {
      stop(
        "Predicted probability matrix is missing columns: ",
        paste(missing, collapse = ", ")
      )
    }
  }

  p[, levels, drop = FALSE]
}
