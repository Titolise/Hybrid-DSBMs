# ==============================================================================
# Simulation Benchmark: Hybrid Dynamic Stochastic Block Model (H-DSBM)
# ==============================================================================

library(hdsbm)
library(mclust)

# ------------------------------------------------------------------------------
# 1. Experimental Configuration
# ------------------------------------------------------------------------------
mod     <- 1          # Scenario: 1 = Assortative (High intra), 2 = Disassortative (High inter)
n       <- 20         # Number of vertices (nodes)
k       <- 3          # Number of latent communities
TT      <- 6          # Number of discrete time snapshots
nrep    <- 15         # Number of perturbed SVD multi-start restarts
persist <- "high"     # Markov chain persistence ("high" vs. "low")
seed_sim <- 42        # Master seed for the simulation run

# ------------------------------------------------------------------------------
# 2. Model Parameters Setup
# ------------------------------------------------------------------------------
# 2.1 Initial class distribution (uniform across latent blocks)
piv0 <- rep(1 / k, k)

# 2.2 First-order Markov transition matrix (Pi)
# High persistence implies lower probability of transitioning between classes
rho <- if (persist == "high") 0.1 else 0.3

# Vectorized Toeplitz-like construction based on community distance
Pi0 <- rho^abs(outer(seq_len(k), seq_len(k), "-"))
Pi0 <- sweep(Pi0, 1, rowSums(Pi0), "/")  # Row-normalize to valid stochastic matrix

# 2.3 Dyadic connection probability matrix (Psi)
set.seed(6)  # Seed for reproducible ground-truth block parameters

intra_val <- if (mod == 1) 0.30 else 0.03
inter_val <- if (mod == 1) 0.03 else 0.30

Psi0       <- matrix(inter_val, nrow = k, ncol = k)
diag(Psi0) <- intra_val * runif(k, min = 0.5, max = 1.5)

# ------------------------------------------------------------------------------
# 3. Dynamic Network Generation
# ------------------------------------------------------------------------------
cat("[Simulation] Generating temporal network under DSBM...\n")
sim_data <- draw_sn_dyn(
  n   = n,
  k   = k,
  TT  = TT,
  piv = piv0,
  Pi  = Pi0,
  Psi = Psi0
)

# Extract ground-truth community trajectories (converted to 1-indexed)
U_true <- sim_data$U + 1

# ------------------------------------------------------------------------------
# 4. Model Estimation via Perturbed SVD Multi-Start
# ------------------------------------------------------------------------------
cat("[Estimation] Fitting H-DSBM with perturbed spectral restarts...\n")
modello <- hdsbm_fit(
  Y      = sim_data$Y,
  k      = k,
  nrep   = nrep,
  maxit  = 2000,
  tol    = 1e-10,
  verbose = TRUE,
  seed   = seed_sim
)

# ------------------------------------------------------------------------------
# 5. Performance Evaluation
# ------------------------------------------------------------------------------
# Evaluate clustering agreement using Adjusted Rand Index (ARI)
ari_greedy  <- mclust::adjustedRandIndex(modello$best_fit$cl,  U_true)
ari_viterbi <- mclust::adjustedRandIndex(modello$best_fit$clv, U_true)

cat("\n====================================================\n")
cat("                BENCHMARK RESULTS                   \n")
cat("====================================================\n")
print(modello)
cat("----------------------------------------------------\n")
cat(sprintf("ARI (Pointwise Greedy): %.4f\n", ari_greedy))
cat(sprintf("ARI (Global Viterbi):   %.4f\n", ari_viterbi))
cat("====================================================\n")