# hdsbm: Hybrid Dynamic Stochastic Block Models

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/gianfilippotito/hdsbm)
[![R Version](https://img.shields.io/badge/R-%3E%3D%203.5-blue.svg)](https://www.r-project.org/)

An R package providing fast estimation, simulation, and model selection routines for discrete-time **Dynamic Stochastic Block Models (DSBM)**. The package features an accelerated **Classification Expectation-Maximization (CEM)** algorithm, an asymmetric M-step, high-performance C++ subroutines powered by `Rcpp` and `Armadillo`, and an innovative **Perturbed SVD Multi-Start Initialization Strategy**.

Developed as part of a Master's Thesis in Statistics by **Gianfilippo Tito**.

---

## 📌 Theoretical Framework

Dynamic Stochastic Block Models extend static stochastic block models to temporal graph structures, allowing the latent block memberships of nodes to transition across discrete time snapshots $t = 1, \dots, T$ governed by a first-order Markov chain.

Let $\mathbf{Y} \in \{0, 1\}^{n \times n \times T}$ denote a binary temporal adjacency tensor observed over $n$ vertices across $T$ discrete observation horizons. Let $U_{i, t} \in \{1, \dots, k\}$ be the latent community membership of node $i$ at time $t$. 

The model parameters are parameterized by $\boldsymbol{\theta} = (\boldsymbol{\pi}_0, \boldsymbol{\Pi}, \boldsymbol{\Psi})$:

* **Initial class proportions**: $\boldsymbol{\pi}_0 \in \Delta^{k-1}$, where $\pi_{0, u} = \mathbb{P}(U_{i, 1} = u)$;
* **First-order Markov transition matrix**: $\boldsymbol{\Pi} \in [0, 1]^{k \times k}$, where $\Pi_{u, v} = \mathbb{P}(U_{i, t} = v \mid U_{i, t-1} = u)$ and $\sum_{v=1}^k \Pi_{u, v} = 1$;
* **Dyadic connection probability matrix**: $\boldsymbol{\Psi} \in [0, 1]^{k \times k}$, where $\Psi_{u, v} = \mathbb{P}(Y_{i, j, t} = 1 \mid U_{i, t} = u, U_{j, t} = v)$.

Due to the combinatorial explosion of latent trajectory paths and the multi-modal topography of the complete-data log-likelihood surface, standard local search heuristics frequently get trapped in sub-optimal local extrema. `hdsbm` addresses this problem by coupling a **Hybrid CEM algorithm** with a structured **Perturbed SVD Spectral Multi-Start** strategy.

---

## 🔬 Perturbed SVD Multi-Start Strategy

A key contribution of the package is the **Perturbed SVD Multi-Start scheme**, combining the asymptotic properties of spectral graph clustering with controlled stochastic exploration of the likelihood landscape.

```
          Temporal Network Tensor Y (n x n x T)
                            │
                            ▼
           Time-Averaged Adjacency Matrix Y_mean
                            │
                            ▼
              Truncated SVD Decomposition
                            │
                            ▼
             K-Means Base Partition (cl_base)
                            │
          ┌─────────────────┴─────────────────┐
          ▼                                   ▼
    [Base Start]               [h = 1, ..., nrep Perturbations]
     Tau ~ K-means             Perturbation intensity α ~ U(0.05, 0.35)
                               Stochastic class reassignment from cl_base
                               Column-wise normalization -> Taur
                                      │
                                      ▼
                             Hybrid CEM Fit with Taur
                                      │
          └─────────────────┬─────────────────┘
                            ▼
              Select Optimal Estimate Based on
                Complete-Data Log-Likelihood
```

### Algorithmic Workflow

1. **Spectral Base Partitioning**:
   * Compute the time-averaged network matrix:
     $$
     \bar{\mathbf{Y}} = \frac{1}{T} \sum_{t=1}^T \mathbf{Y}_t
     $$
   * Extract the Singular Value Decomposition (SVD) of $\bar{\mathbf{Y}}$ and retain the top $k$ left singular vectors $\mathbf{U}_k \in \mathbb{R}^{n \times k}$:
     $$
     \bar{\mathbf{Y}} \approx \mathbf{U}_k \boldsymbol{\Sigma}_k \mathbf{V}_k^\top
     $$
   * Run $k$-means on the rows of $\mathbf{U}_k$ to construct the deterministic baseline partition $\mathbf{c}^{\text{base}} \in \{1, \dots, k\}^n$.

2. **Controlled Stochastic Perturbations**:
   * For each restart $h \in \{1, \dots, \text{nrep}\}$:
     * Sample a random perturbation fraction $\alpha_h \sim \mathcal{U}(0.05, 0.35)$.
     * For each node $i$ and snapshot $t$, retain $c_i^{\text{base}}$ with probability $1 - \alpha_h$, or reassign it to a uniformly chosen alternative block $v \in \{1, \dots, k\} \setminus \{c_i^{\text{base}}\}$ with probability $\alpha_h$.
     * Construct a regularized posterior tensor $\boldsymbol{\tau}^{(h)} \in [0, 1]^{k \times n \times T}$ ($0.98$ mass on the selected cluster, $0.01$ distributed uniformly across remainder), followed by column normalization:
       $$
       \sum_{u=1}^k \tau_{u, i, t}^{(h)} = 1
       $$

3. **Global Model Selection**:
   * Execute the Hybrid CEM engine starting from the baseline and each perturbed state $\boldsymbol{\tau}^{(h)}$.
   * Retain the solution maximizing the complete-data log-likelihood:
     $$
     \hat{\boldsymbol{\theta}} = \arg\max_{h \in \{0, 1, \dots, \text{nrep}\}} \ell\left(\boldsymbol{\theta}^{(h)} \mid \mathbf{Y}\right)
     $$

---

## ⚡ Core Engine Highlights

* **Forward-Backward Recursions (E-step)**: Computes soft responsibilities over hidden temporal transitions via numerically stabilized Hidden Markov Model forward-backward routines.
* **Asymmetric M-step**: Simultaneously integrates soft continuous responsibilities and hard categorical assignments to update $\boldsymbol{\pi}_0$, $\boldsymbol{\Pi}$, and $\boldsymbol{\Psi}$.
* **C++ Accelerated Combinatorial Greedy Search**: Accelerates dyadic reallocation passes through compiled C++ subroutines (`src/Core_hyb.cpp`) using `RcppArmadillo`.
* **Dual-Path Decoding**: Delivers both pointwise greedy classifications (`cl`) and globally optimal Viterbi state trajectories (`clv`).
* **Integrated Completed Likelihood (ICL)**: Tailored penalization criterion separating structural connectivity parameters from Markov chain dynamic parameters.

---

## 📂 Repository Structure

The package follows standard CRAN directory conventions:

```text
hdsbm/
├── DESCRIPTION             # Package metadata, dependencies, and licensing
├── NAMESPACE               # Exported functions and C++ symbol registration
├── LICENSE / LICENSE.md    # MIT License terms
├── .Rbuildignore           # Rules for packaging and build exclusion
├── .gitignore              # Git ignore configuration for Rcpp binaries and caches
├── data/
│   └── toy_hdsbm.rda       # Built-in synthetic sample dataset
├── data-raw/
│   └── toy_hdsbm.R         # Reproducible data-generation script
├── man/                    # Roxygen2 Rd documentation files
├── R/
│   ├── best_perm.R         # Hungarian-style label alignment for ARI evaluation
│   ├── data.R              # Dataset documentation
│   ├── draw_sn_dyn.R       # Dynamic network generative simulator
│   ├── est_hyb_sbm_dyn_dec.R # CEM algorithm with forward-backward & Viterbi passes
│   ├── hdsbm_fit.R         # Top-level API with SVD multi-start and S3 print methods
│   ├── hdsbm-package.R     # Package-level documentation and Rcpp registration
│   ├── icl_dyn_hyb.R       # Integrated Completed Likelihood calculation
│   └── RcppExports.R       # Automatically generated Rcpp wrappers
├── src/
│   ├── Core_hyb.cpp        # High-performance C++ likelihood & greedy search functions
│   └── RcppExports.cpp     # Compiled Rcpp registration routines
└── Simulation.R            # Benchmarking and Monte Carlo simulation script
```

---

## ⚙️ Requirements & Installation

### Prerequisites

* **R (>= 3.5.0)** (R >= 4.0 recommended).
* A functional C++ compiler with C++11 support:
  * **macOS**: Xcode Command Line Tools (`xcode-select --install`).
  * **Windows**: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version.
  * **Linux**: `g++` or `clang++` (`build-essential` on Debian/Ubuntu).

### Install via `remotes` / `devtools`

You can install `hdsbm` directly from GitHub:

```r
# If not already installed:
install.packages("remotes")

# Install hdsbm directly from GitHub:
remotes::install_github("gianfilippotito/hdsbm")
```

Or clone the repository and build locally:

```bash
git clone https://github.com/gianfilippotito/hdsbm.git
cd hdsbm
R CMD INSTALL .
```

---

## 🚀 Quick Start Guide

### 1. Using the Included Dataset (`toy_hdsbm`)

The package includes a ready-to-use temporal network dataset `toy_hdsbm` ($n = 25$ nodes, $k = 3$ latent classes, $T = 4$ snapshots):

```r
library(hdsbm)
library(mclust)

# 1. Load the dataset (available lazily)
data(toy_hdsbm)

# 2. Fit the model with 5 perturbed spectral restarts
modello <- hdsbm_fit(
  Y     = toy_hdsbm$Y,
  k     = 3,
  nrep  = 5,
  seed  = 123
)

# 3. Print summary results (S3 method)
print(modello)

# 4. Evaluate clustering accuracy vs. ground truth
ari_greedy  <- adjustedRandIndex(modello$best_fit$cl,  toy_hdsbm$U_true)
ari_viterbi <- adjustedRandIndex(modello$best_fit$clv, toy_hdsbm$U_true)

cat(sprintf("Adjusted Rand Index (Pointwise Greedy): %.4f\n", ari_greedy))
cat(sprintf("Adjusted Rand Index (Global Viterbi):   %.4f\n", ari_viterbi))
```

### 2. Simulating and Estimating a Dynamic Network

Generate synthetic dynamic networks and test different persistent Markov regimes:

```r
library(hdsbm)

# Dimensional parameters
n  <- 30     # Nodes
k  <- 3      # Communities
TT <- 5      # Time windows

# Ground-truth parameters
piv0 <- rep(1 / k, k)

# High persistence Markov transition matrix
rho <- 0.1
Pi0 <- rho^abs(outer(seq_len(k), seq_len(k), "-"))
Pi0 <- sweep(Pi0, 1, rowSums(Pi0), "/")

# Assortative connectivity matrix (high within-group, low between-group)
Psi0 <- matrix(0.04, nrow = k, ncol = k)
diag(Psi0) <- c(0.35, 0.40, 0.45)

# Simulate dynamic network
sim <- draw_sn_dyn(n = n, k = k, TT = TT, piv = piv0, Pi = Pi0, Psi = Psi0)

# Estimate with SVD multi-start
fit <- hdsbm_fit(
  Y      = sim$Y,
  k      = k,
  nrep   = 10,
  maxit  = 150,
  seed   = 42,
  verbose = TRUE
)

# Access estimated parameters
print(fit)
cat("Optimal ICL Score:", fit$icl$ICL, "\n")
```

---

## 📖 Function Reference

| Function | Description |
| :--- | :--- |
| `hdsbm_fit(Y, k, nrep, ...)` | **Main API function**. Runs the Hybrid CEM engine using the Perturbed SVD multi-start strategy and returns an S3 `hdsbm` object. |
| `draw_sn_dyn(n, k, TT, piv, Pi, Psi)` | Generates synthetic dynamic network tensors and true latent paths under the DSBM framework. |
| `icl_dyn_hyb(out_pred, n, TT, k)` | Computes the Integrated Completed Likelihood (ICL) with dedicated network and Markov transition penalties. |
| `best_perm(Utrue, Uprop)` | Finds the optimal permutation of class labels to match estimated clusters to ground truth. |
| `est_pred_sbm_dyn_decoding(Y, k, ...)` | Low-level estimation routine executing the Hybrid CEM algorithm with Viterbi decoding. |

---

## 📊 Evaluation & Model Selection

* **Adjusted Rand Index (ARI)**: Evaluates partition recovery accuracy against ground-truth community trajectories, invariant to label switching.
* **Integrated Completed Likelihood (ICL)**:
  $$
  \text{ICL} = \ell_c(\hat{\boldsymbol{\theta}} \mid \mathbf{Y}, \hat{\mathbf{U}}) - \frac{k(k+1)}{4}\log\left(T \frac{n(n-1)}{2}\right) - \frac{k(k-1)}{2}\log\left(n(T-1)\right) - \frac{k-1}{2}\log(n)
  $$
  Penalizes connection probabilities, transition matrices, and initial distributions separately to guarantee consistent model selection.

---

## 📝 License & Author

* **Author**: Gianfilippo Tito ([gianfilippotito@gmail.com](mailto:gianfilippotito@gmail.com))
* **License**: Released under the [MIT License](LICENSE.md).