#' Hybrid DSBM Estimation with Perturbed SVD Multi-Start
#'
#' Estimates a Dynamic Stochastic Block Model (DSBM) by combining
#' the Hybrid CEM algorithm with a multi-start strategy based on the spectral
#' decomposition (SVD) of the time-averaged adjacency matrix and stochastic perturbations.
#'
#' @param Y Binary 3D array of dimensions \eqn{n \times n \times T} (the temporal adjacency matrix).
#' @param k Integer, number of latent communities/blocks.
#' @param nrep Integer, number of random restarts perturbed via SVD (default: 10). If 0, only the initial fit is performed.
#' @param perturb_range Numeric vector of length 2 indicating the lower and upper limits
#'   of the random perturbation rate (default: \code{c(0.05, 0.35)}).
#' @param maxit Integer, maximum number of CEM iterations per restart (default: 200).
#' @param tol Numeric, convergence tolerance for the log-likelihood (default: 1e-10).
#' @param verbose Logical; if \code{TRUE}, prints restart progress to the console (default: \code{TRUE}).
#' @param seed Optional integer to ensure reproducibility.
#'
#' @return An object of class \code{hdsbm} containing:
#' \itemize{
#'   \item \code{best_fit}: the estimate with the highest log-likelihood (parameters \eqn{\pi_0, \Pi, \Psi}, classifications \code{cl} and \code{clv}).
#'   \item \code{lk_trace}: vector containing the log-likelihoods obtained from each restart. #'   \item \code{best_restart}: the index of the restart that achieved the optimal value (0 = baseline).
#'   \item \code{icl}: ICL score calculated on the optimal fit.
#'   \item \code{params}: model dimensional parameters (\code{n}, \code{k}, \code{TT}).
#' }
#'
#' @importFrom stats kmeans runif
#' @export
 
hdsbm_fit <- function(Y,
                      k,
                      nrep = 10,
                      perturb_range = c(0.05, 0.35),
                      maxit = 200,
                      tol = 1e-10,
                      verbose = TRUE,
                      seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Controlli dimensionali e validazione input
  if (!is.array(Y) || length(dim(Y)) != 3) {
    stop("`Y` must be a three-dimensional array (n x n x TT).")
  }
  n  <- dim(Y)[1]
  TT <- dim(Y)[3]
  if (dim(Y)[2] != n) {
    stop("The time slices of `Y` must be square matrices (n x n).")
  }
  if (k <= 0 || k > n) {
    stop("`k` It must be a positive integer between 1 and n.")
  }
  
  # Gestione self-loop (impostati a NA come nello script di simulazione)
  YY <- Y
  for (t in seq_len(TT)) {
    diag(YY[, , t]) <- NA
  }
  
  # 2. Fit Iniziale Baseline (K-means sulle fette concatenate)
  if (verbose) cat("[hdsbm] Initial baseline fit startup...\n")
  best_fit <- est_pred_sbm_dyn_decoding(
    YY,
    k     = k,
    start = 0,
    tol   = tol,
    maxit = maxit
  )
  
  best_lk      <- best_fit$lk
  lk_trace     <- c(baseline = best_lk)
  best_restart <- 0
  
  # 3. Multi-Start SVD Perturbato
  if (nrep > 0) {
    if (verbose) {
      cat(sprintf("[hdsbm] Execution of perturbed %d restarts (SVD-based)...\n", nrep))
    }
    
    # Partizione spettrale di base dalla matrice media nel tempo
    Y_mean <- apply(YY, c(1, 2), mean, na.rm = TRUE)
    Y_mean[is.na(Y_mean)] <- 0
    
    # Truncated SVD sui primi k autovettori sinistri
    svd_decomp <- svd(Y_mean)
    vecs <- svd_decomp$u[, seq_len(k), drop = FALSE]
    cl_base <- stats::kmeans(vecs, centers = k, nstart = 50)$cluster
    
    # Loop sui restart
    for (h in seq_len(nrep)) {
      perturb_rate <- stats::runif(1, min = perturb_range[1], max = perturb_range[2])
      Taur <- array(0.01, c(k, n, TT))
      
      for (t in seq_len(TT)) {
        for (i in seq_len(n)) {
          assigned_class <- cl_base[i]
          if (stats::runif(1) < perturb_rate) {
            alt_class <- setdiff(seq_len(k), cl_base[i])
            if (length(alt_class) > 0) {
              assigned_class <- sample(alt_class, 1)
            }
          }
          Taur[assigned_class, i, t] <- 0.98
        }
        # Normalizzazione stocastica delle colonne
        Taur[, , t] <- sweep(Taur[, , t], 2, colSums(Taur[, , t]), "/")
      }
      
      if (verbose) {
        cat(sprintf("  Restart %d/%d (Perturb: %.1f%%)... ", h, nrep, perturb_rate * 100))
      }
      
      fit_h <- est_pred_sbm_dyn_decoding(
        YY,
        k     = k,
        start = 2,
        tol   = tol,
        Tau   = Taur,
        maxit = maxit
      )
      
      lk_trace[paste0("restart_", h)] <- fit_h$lk
      
      if (fit_h$lk > best_lk) {
        if (verbose) cat(sprintf("-> NEW OPTIMAL! (lk = %.3f)\n", fit_h$lk))
        best_lk      <- fit_h$lk
        best_fit     <- fit_h
        best_restart <- h
      } else {
        if (verbose) cat(sprintf("(lk = %.3f)\n", fit_h$lk))
      }
    }
  }
  
  # 4. Calcolo metrica ICL finale
  icl_eval <- icl_dyn_hyb(best_fit, n = n, TT = TT, k = k)
  
  # 5. Output strutturato con classe S3
  res <- list(
    best_fit     = best_fit,
    lk_trace     = lk_trace,
    best_restart = best_restart,
    icl          = icl_eval,
    params       = list(n = n, k = k, TT = TT, nrep = nrep)
  )
  
  class(res) <- "hdsbm"
  return(res)
}


#' @export
print.hdsbm <- function(x, ...) {
  cat("====================================================\n")
  cat("  Hybrid Dynamic Stochastic Block Model (H-DSBM)\n")
  cat("====================================================\n")
  cat(sprintf("Nodes (n): %d  |  Communities (k): %d  |  Snapshots (T): %d\n",
              x$params$n, x$params$k, x$params$TT))
  cat(sprintf("Total Restart: %d  |  Winning Restart: %s\n",
              x$params$nrep,
              if (x$best_restart == 0) "Baseline (K-means)" else paste0("Restart #", x$best_restart)))
  cat(sprintf("Optimal Log-Likelihood: %.4f\n", x$best_fit$lk))
  cat(sprintf("ICL Score: %.2f\n", x$icl$ICL))
  cat("----------------------------------------------------\n")
  cat("Estimated connection matrix (Psi):\n")
  print(round(x$best_fit$Psi, 4))
  cat("\nEstimated transition matrix (Pi):\n")
  print(round(x$best_fit$Pi, 4))
  invisible(x)
}
