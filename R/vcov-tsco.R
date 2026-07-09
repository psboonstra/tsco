#' Variance-covariance matrix for a TSCO model
#'
#' @param object An object of class `"tsco"`.
#' @param stage Which stage to extract: `"both"`, `"stage1"`, or `"stage2"`.
#' @param ... Ignored.
#'
#' @return A variance-covariance matrix. If `stage = "both"`, the returned
#'   matrix is block diagonal with stage-specific coefficient names prefixed by
#'   `"stage1:"` and `"stage2:"`.
#'
#' @method vcov tsco
#' @export
vcov.tsco <- function(object, stage = c("both", "stage1", "stage2"), ...) {
  stage <- match.arg(stage)

  V1 <- as.matrix(stats::vcov(object$fit_stage1))
  V2 <- as.matrix(stats::vcov(object$fit_stage2))

  if (stage == "stage1") {
    return(V1)
  }

  if (stage == "stage2") {
    return(V2)
  }

  p1 <- nrow(V1)
  p2 <- nrow(V2)

  if (is.null(p1) || is.null(p2) || p1 < 1L || p2 < 1L) {
    stop("Could not extract valid stage-specific covariance matrices.",
         call. = FALSE)
  }

  i1 <- seq_len(p1)
  i2 <- p1 + seq_len(p2)

  V <- matrix(0, nrow = p1 + p2, ncol = p1 + p2)

  V[i1, i1] <- V1
  V[i2, i2] <- V2

  rn1 <- rownames(V1)
  rn2 <- rownames(V2)

  if (is.null(rn1)) {
    rn1 <- names(stats::coef(object$fit_stage1))
  }

  if (is.null(rn2)) {
    rn2 <- names(stats::coef(object$fit_stage2))
  }

  rn <- c(
    paste0("stage1:", rn1),
    paste0("stage2:", rn2)
  )

  dimnames(V) <- list(rn, rn)

  V
}
