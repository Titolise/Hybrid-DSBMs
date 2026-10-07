icl_dyn_hyb <- function(out_pred, n, TT, k) {
  
  # out_pred: output returned by est_pred_sbm_dyn_decoding()
  # n: number of nodes
  # TT: number of time steps
  # k: number of latent classes (Q)
  
  lk <- out_pred$lk
  
  # Psi Connectivity Matrix Penalty
  n_dyads_total <- TT * (n * (n - 1) / 2)
  pen_psi <- 0.5 * (k * (k + 1) / 2) * log(n_dyads_total)
  
  # Rho Transition Matrix Penalty
  n_transitions_total <- n * (TT - 1)
  pen_rho <- 0.5 * (k * (k - 1)) * log(n_transitions_total)
  
  # Penalty for initial PIV proportions
  pen_pi <- 0.5 * (k - 1) * log(n)
  
  # ICL
  icl <- lk - (pen_psi + pen_rho + pen_pi)
  
  return(list(
    ICL = icl,
    lk = lk,
    pen_total = pen_psi + pen_rho + pen_pi,
    pen_components = c(psi = pen_psi, rho = pen_rho, pi = pen_pi)
  ))
}