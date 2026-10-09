#' Integrated Completed Likelihood for Hybrid DSBM
#'
#' Calculates the ICL model selection criterion by penalizing both
#' connection probabilities and Markovian transition probabilities.
#'
#' @param out_pred List object returned by \code{est_pred_sbm_dyn_decoding}.
#' @param n Integer, number of nodes.
#' @param TT Integer, number of time points.
#' @param k Integer, number of latent communities.
#'
#' @return List containing the \code{ICL} value, the log-likelihood \code{lk},
#' and the penalty decomposition.
#' @export

icl_dyn_hyb <- function(out_pred, n, TT, k) {
  lk <- out_pred$lk
  
  # Psi connection matrix penalty
  n_dyads_total <- TT * (n * (n - 1) / 2)
  pen_psi <- 0.5 * (k * (k + 1) / 2) * log(n_dyads_total)
  
  # Pi transition matrix penalty
  n_transitions_total <- n * (TT - 1)
  pen_rho <- 0.5 * (k * (k - 1)) * log(n_transitions_total)
  
  # Initial proportions piv penalty
  pen_pi <- 0.5 * (k - 1) * log(n)
  
  pen_tot <- pen_psi + pen_rho + pen_pi
  icl <- lk - pen_tot
  
  list(
    ICL = icl,
    lk = lk,
    pen_total = pen_tot,
    pen_components = c(psi = pen_psi, rho = pen_rho, pi = pen_pi)
  )
}