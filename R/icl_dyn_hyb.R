#' Integrated Completed Likelihood (ICL) for Dynamic SBM
#'
#' Computes the Integrated Completed Likelihood (ICL) criterion for
#' Dynamic Stochastic Block Models. Penalizes dyadic connection probabilities,
#' first-order Markov transition parameters, and initial community proportions.
#'
#' @param object An object of class \code{hdsbm} (returned by \code{\link{hdsbm_fit}}),
#'   or a list returned by \code{\link{est_pred_sbm_dyn_decoding}}.
#' @param n Optional integer, number of nodes. Inferred automatically if \code{NULL}.
#' @param TT Optional integer, number of time points. Inferred automatically if \code{NULL}.
#' @param k Optional integer, number of latent classes. Inferred automatically if \code{NULL}.
#' @param detailed Logical; if \code{FALSE} (default), returns a single numeric value (like \code{AIC} or \code{BIC}).
#'   If \code{TRUE}, returns a list with detailed penalty decompositions.
#' @param ... Additional arguments passed to methods.
#'
#' @return If \code{detailed = FALSE}, a numeric value representing the ICL score.
#'   If \code{detailed = TRUE}, a list containing:
#'   \itemize{
#'     \item \code{ICL}: The total ICL score.
#'     \item \code{lk}: The complete-data log-likelihood.
#'     \item \code{pen_total}: The sum of all penalties.
#'     \item \code{pen_components}: Vector with individual penalties (\code{psi}, \code{rho}, \code{pi}).
#'   }
#'
#' @export
#' @examples
#' data(toy_hdsbm)
#' model <- hdsbm_fit(toy_hdsbm$Y, k = 3, nrep = 2, maxit = 20, verbose = FALSE)
#'
#' # Standard numeric value (ideal for model comparison)
#' icl(model)
#'
#' # Detailed penalty breakdown
#' icl(model, detailed = TRUE)
icl <- function(object, n = NULL, TT = NULL, k = NULL, detailed = FALSE, ...) {
  UseMethod("icl")
}

#' @rdname icl
#' @export
icl.hdsbm <- function(object, n = NULL, TT = NULL, k = NULL, detailed = FALSE, ...) {
  n_val  <- n  %||% object$params$n
  tt_val <- TT %||% object$params$TT
  k_val  <- k  %||% object$params$k
  lk_val <- object$best_fit$lk
  
  res <- compute_icl_values(lk = lk_val, n = n_val, TT = tt_val, k = k_val)
  
  if (detailed) res else res$ICL
}

#' @rdname icl
#' @export
icl.default <- function(object, n = NULL, TT = NULL, k = NULL, detailed = FALSE, ...) {
  if (!is.list(object) || is.null(object$lk)) {
    stop("`object` must be an `hdsbm` object or a list containing `$lk` from `est_pred_sbm_dyn_decoding()`.")
  }
  
  n_val  <- n  %||% (if (!is.null(object$cl)) nrow(object$cl) else NULL)
  tt_val <- TT %||% (if (!is.null(object$cl)) ncol(object$cl) else NULL)
  k_val  <- k  %||% (if (!is.null(object$piv)) length(object$piv) else if (!is.null(object$Psi)) nrow(object$Psi) else NULL)
  
  if (is.null(n_val) || is.null(tt_val) || is.null(k_val)) {
    stop("Could not infer all dimensions (n, TT, k) from `object`. Please supply them explicitly.")
  }
  
  res <- compute_icl_values(lk = object$lk, n = n_val, TT = tt_val, k = k_val)
  
  if (detailed) res else res$ICL
}

#' @rdname icl
#' @export
ICL <- function(object, ...) {
  icl(object, ...)
}

#' @rdname icl
#' @param out_pred Alias for \code{object} for backward compatibility.
#' @export
icl_dyn_hyb <- function(out_pred, n = NULL, TT = NULL, k = NULL) {
  icl(object = out_pred, n = n, TT = TT, k = k, detailed = TRUE)
}

# --- Internal Helper Routine ---
compute_icl_values <- function(lk, n, TT, k) {
  # Psi connection matrix penalty
  n_dyads_total <- TT * (n * (n - 1) / 2)
  pen_psi <- 0.5 * (k * (k + 1) / 2) * log(n_dyads_total)
  
  # Pi transition matrix penalty
  n_transitions_total <- n * (TT - 1)
  pen_rho <- 0.5 * (k * (k - 1)) * log(n_transitions_total)
  
  # Initial proportions piv penalty
  pen_pi <- 0.5 * (k - 1) * log(n)
  
  pen_tot <- pen_psi + pen_rho + pen_pi
  icl_val <- lk - pen_tot
  
  list(
    ICL            = icl_val,
    lk             = lk,
    pen_total      = pen_tot,
    pen_components = c(psi = pen_psi, rho = pen_rho, pi = pen_pi)
  )
}

# Null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b