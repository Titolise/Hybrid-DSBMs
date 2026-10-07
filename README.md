# Dynamic Stochastic Block Model (DSBM) – Hybrid CEM Estimation

This repository contains the **R** and **C++** implementation (via `Rcpp` / `RcppArmadillo`) of the **Hybrid Classification Expectation-Maximization (CEM)** estimation algorithm for discrete-time **Dynamic Stochastic Block Models (DSBM)**, developed as part of a Master's Thesis in Statistics.

## 📌 Model Overview

The model analyzes sequences of binary, undirected networks represented as adjacency tensors $\mathbf{Y} \in \{0, 1\}^{n \times n \times T}$, where:
* $n$ is the number of nodes;
* $T$ is the number of time horizons (slices);
* $k$ is the number of latent communities (blocks).

Node community membership evolves over time according to a **first-order Markov chain**:
* $\boldsymbol{\pi}_0 \in \Delta^{k-1}$: initial class probabilities at time $t = 1$;
* $\boldsymbol{\Pi} \in \mathbb{R}^{k \times k}$: first-order Markov transition probability matrix;
* $\boldsymbol{\Psi} \in [0, 1]^{k \times k}$: symmetric edge density matrix (intra- and inter-community connection probabilities).

### Hybrid CEM Algorithm

Unlike standard Variational EM (VEM) or naive CEM formulations:
1. **E-step**: Evaluates posterior probabilities over node trajectories using forward-backward recursive formulations (Hidden Markov Model framework).
2. **M-step**: Analytically updates transition parameters $\boldsymbol{\pi}_0$, $\boldsymbol{\Pi}$, and the dyadic connectivity matrix $\boldsymbol{\Psi}$.
3. **Combinatorial Optimization (C++)**: Performs a local greedy search over discrete cluster allocations to maximize the conditional complete-data log-likelihood, accelerated via `RcppArmadillo`.
4. **Global Decoding (Viterbi)**: Computes the optimal sequential trajectory of cluster assignments using the Viterbi algorithm (`clv`), alongside the pointwise greedy classification output (`cl`).

## 📂 Repository Structure

```
.
├── Cpp/
│   └── Core_hyb.cpp           # High-performance C++ routines (complete log-likelihood & greedy search)
├── src/
│   ├── best_perm.R            # Hungarian/permutation alignment to resolve label switching
│   ├── draw_sn_dyn.R          # Synthetic dynamic network generator under DSBM
│   ├── est_hyb_sbm_dyn_dec.R  # Main Hybrid DSBM estimation procedure with Viterbi decoding
│   └── icl_dyn_hyb.R          # Model selection criterion (Integrated Completed Likelihood - ICL)
├── Simulation.R               # Monte Carlo benchmark script with random spectral perturbations
└── README.md                  # Project documentation
```

## ⚙️ Requirements & Dependencies

The codebase requires **R** (>= 4.0.0) and a C++11 compliant compiler (e.g., `g++`, `clang`, or Rtools for Windows users).

### Required R Packages:
```r
install.packages(c(
  "Rcpp",
  "RcppArmadillo",
  "mclust",
  "gtools"
))
```

## 🚀 Usage

### 1. Running the Simulation Study

To execute the complete Monte Carlo benchmark pipeline:
```r
source("Simulation.R")
```

Key simulation parameters can be configured directly inside `Simulation.R`:
* `n`: Number of nodes in the network (default: `20`);
* `TT`: Number of discrete time steps (default: `6`);
* `k`: Number of latent classes (default: `3`);
* `B`: Number of Monte Carlo replicates (default: `100`);
* `nrep`: Number of perturbed spectral restarts to mitigate local maxima;
* `persist`: Markov chain persistence degree (`"high"` vs. `"low"`).

### 2. Minimal Working Example (Synthetic Data Estimation)

```r
library(Rcpp)
library(RcppArmadillo)
library(mclust)

# Load routines
sourceCpp("./Cpp/Core_hyb.cpp")
source("./src/draw_sn_dyn.R")
source("./src/est_hyb_sbm_dyn_dec.R")
source("./src/icl_dyn_hyb.R")

# 1. Simulate a 3-class DSBM with 30 nodes over 5 time steps
n <- 30
k <- 3
TT <- 5
piv0 <- rep(1 / k, k)
Pi0 <- diag(0.8, k) + matrix(0.2 / (k - 1), k, k)
diag(Pi0) <- 0.8
Psi0 <- matrix(0.05, k, k)
diag(Psi0) <- 0.35

sim <- draw_sn_dyn(n, k, TT, piv0, Pi0, Psi0)
Y <- sim$Y

# 2. Set diagonal entries to NA (no self-loops)
for (t in 1:TT) diag(Y[,, t]) <- NA

# 3. Fit the Hybrid DSBM
fit <- est_pred_sbm_dyn_decoding(Y, k = k, start = 0, maxit = 200)

# 4. Evaluate clustering accuracy via Adjusted Rand Index (ARI)
ari_greedy  <- adjustedRandIndex(fit$cl,  sim$U + 1)
ari_viterbi <- adjustedRandIndex(fit$clv, sim$U + 1)

cat(sprintf("ARI (Greedy): %.3f | ARI (Viterbi): %.3f\n", ari_greedy, ari_viterbi))

# 5. Model selection via dynamic ICL
icl_res <- icl_dyn_hyb(fit, n = n, TT = TT, k = k)
cat(sprintf("Dynamic ICL value: %.2f\n", icl_res$ICL))
```

## 📊 Evaluation Metrics

* **Adjusted Rand Index (ARI)**: Assesses clustering concordance with true latent state trajectories $U_{\text{true}}$, inherently invariant to label switching.
* **Best Permutation Alignment**: Finds the optimal bijective mapping across cluster labels to compute contingency tables and classification accuracy rates.
* **Integrated Completed Likelihood (ICL)**: Tailored for dynamic networks with penalizations separating Markov transition parameters ($\boldsymbol{\Pi}$) and dyadic edge probabilities ($\boldsymbol{\Psi}$).
