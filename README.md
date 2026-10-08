# Hybrid Dynamic Stochastic Block Models (H-DSBM)

Implementation and Monte Carlo benchmarking framework for discrete-time Dynamic Stochastic Block Models via an accelerated **Classification Expectation-Maximization (CEM)** algorithm featuring an asymmetric M-step, high-performance C++ subroutines, and a **Perturbed SVD Multi-Start Initialization Strategy**.

Developed as part of a Master's Thesis in Statistics.

## 📌 Theoretical Framework

Dynamic Stochastic Block Models (DSBM) extend classic stochastic block models to temporal networks, allowing the latent community memberships of nodes to evolve across discrete time slices $t = 1, \dots, T$ according to a first-order Markov chain.

Let $\mathbf{Y} \in \{0, 1\}^{n \times n \times T}$ be the binary temporal adjacency tensor observed on $n$ vertices across $T$ time horizons. Let $U_{i, t} \in \{1, \dots, k\}$ represent the latent community label of node $i$ at time $t$. The model parameters comprise:

* **Initial class proportions**: $\boldsymbol{\pi}_0 \in \Delta^{k-1}$, where $\pi_{0, u} = \mathbb{P}(U_{i, 1} = u)$;

* **First-order Markov transition matrix**: $\boldsymbol{\Pi} \in [0, 1]^{k \times k}$, with $\Pi_{u, v} = \mathbb{P}(U_{i, t} = v \mid U_{i, t-1} = u)$;

* **Connection probability matrix**: $\boldsymbol{\Psi} \in [0, 1]^{k \times k}$, where $\Psi_{u, v} = \mathbb{P}(Y_{i, j, t} = 1 \mid U_{i, t} = u, U_{j, t} = v)$.

Due to the combinatorial space of latent state trajectories and non-convex likelihood surfaces, traditional iterative estimators are prone to getting trapped in suboptimal local modes. This implementation couples a **Hybrid CEM algorithm** with a structured **Perturbed Spectral Multi-Start** strategy to ensure robust global convergence.

## 🔬 Perturbed SVD Multi-Start Strategy

A key contribution of this project is the **Perturbed SVD Multi-Start** scheme. It combines the asymptotic consistency of spectral graph clustering with controlled stochastic exploration of the parameter landscape.

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

### Algorithmic Workflow:

1. **Spectral Base Partitioning**:

   * Compute the time-averaged adjacency matrix across all time frames:
     

     $$
     \bar{\mathbf{Y}} = \frac{1}{T} \sum_{t=1}^T \mathbf{Y}_t
     $$

   * Compute the Singular Value Decomposition (SVD) of $\bar{\mathbf{Y}}$ and retain the top $k$ left singular vectors $\mathbf{U}_k \in \mathbb{R}^{n \times k}$:
     

     $$
     \bar{\mathbf{Y}} \approx \mathbf{U}_k \boldsymbol{\Sigma}_k \mathbf{V}_k^\top
     $$

   * Run $k$-means on the rows of $\mathbf{U}_k$ to obtain a deterministic baseline partition $\mathbf{c}^{\text{base}} \in \{1, \dots, k\}^n$.

2. **Controlled Stochastic Perturbations**:

   * For each restart $h \in \{1, \dots, \text{nrep}\}$:

     * Draw a random perturbation rate $\alpha_h \sim \mathcal{U}(0.05, 0.35)$.

     * For each node $i$ and time slice $t$, retain its baseline assignment $c_i^{\text{base}}$ with probability $1 - \alpha_h$, or swap it to a uniformly drawn alternate cluster $v \in \{1, \dots, k\} \setminus \{c_i^{\text{base}}\}$ with probability $\alpha_h$.

     * Construct a regularized initial posterior tensor $\boldsymbol{\tau}^{(h)} \in [0, 1]^{k \times n \times T}$ (placing $0.98$ mass on the selected cluster, $0.01$ baseline on remaining classes), followed by column normalization:
       

       $$
       \sum_{u=1}^k \tau_{u, i, t}^{(h)} = 1
       $$

3. **Likelihood Maximization Selection**:

   * The Hybrid CEM estimator is executed from both the baseline partition and each perturbed state $\boldsymbol{\tau}^{(h)}$.

   * The final output corresponds to the restart achieving the maximum complete-data log-likelihood:
     

     $$
     \hat{\boldsymbol{\theta}} = \arg\max_{h \in \{0, 1, \dots, \text{nrep}\}} \ell\left(\boldsymbol{\theta}^{(h)} \mid \mathbf{Y}\right)
     $$

## ⚡ Hybrid CEM Algorithm Details

The estimation pipeline combines:

* **Forward-Backward Recursions (E-step)**: Computes node responsibilities over temporal state transitions using Hidden Markov Model forward-backward passes.

* **Asymmetric M-step**: Jointly leverages soft posterior probabilities and hard cluster assignments to update transition parameters ($\boldsymbol{\pi}_0, \boldsymbol{\Pi}$) and the dyadic connectivity matrix $\boldsymbol{\Psi}$.

* **Greedy Combinatorial Search (C++)**: Accelerates node reallocations via an optimized C++ implementation using `RcppArmadillo` (`Cpp/Core_hyb.cpp`).

* **Viterbi Decoding**: Computes the optimal global sequence of hidden states across time (`clv`), alongside pointwise greedy classifications (`cl`).

* **Dynamic ICL Selection**: Computes the Integrated Completed Likelihood (`src/icl_dyn_hyb.R`) with dedicated penalties separating Markov transitions and network edge observations.

## 📂 Repository Structure

```
.
├── Cpp/
│   └── Core_hyb.cpp           # High-performance C++ routines (complk, greedy search)
├── src/
│   ├── best_perm.R            # Hungarian-style label alignment for ARI calculations
│   ├── draw_sn_dyn.R          # Synthetic dynamic network generator under DSBM
│   ├── est_hyb_sbm_dyn_dec.R  # Hybrid CEM estimation routine with Viterbi decoding
│   └── icl_dyn_hyb.R          # Model selection criterion (Integrated Completed Likelihood - ICL)
├── Simulation.R               # Monte Carlo benchmark script with perturbed SVD restarts
├── .gitignore                 # Ignores system files (.DS_Store), compiled binaries, and .RData
└── README.md                  # Project documentation

```

## ⚙️ Requirements & Installation

The project requires **R (>= 4.0.0)** and a C++11 compliant compiler (`g++`, `clang`, or Rtools on Windows).

Install the required R dependencies:

```
install.packages(c(
  "mclust",
  "gtools",
  "Rcpp",
  "RcppArmadillo"
))

```

## 🚀 Usage

### 1. Minimal Working Example (Single Fit)

```
library(Rcpp)
library(RcppArmadillo)
library(mclust)
library(gtools)

# Load core modules
sourceCpp("./Cpp/Core_hyb.cpp")
source("./src/draw_sn_dyn.R")
source("./src/est_hyb_sbm_dyn_dec.R")
source("./src/icl_dyn_hyb.R")

# 1. Setup simulation parameters
n  <- 30     # Number of nodes
k  <- 3      # Number of latent communities
TT <- 5      # Number of time periods

piv0 <- rep(1 / k, k)
Pi0  <- matrix(0.1, k, k)
diag(Pi0) <- 0.8
Pi0  <- diag(1 / rowSums(Pi0)) %*% Pi0

Psi0 <- matrix(0.05, k, k)
diag(Psi0) <- 0.40

# 2. Simulate synthetic temporal network
sim_data <- draw_sn_dyn(n, k, TT, piv0, Pi0, Psi0)
Y <- sim_data$Y

# Set self-loops to NA
for (t in 1:TT) diag(Y[,, t]) <- NA

# 3. Fit Hybrid DSBM
fit <- est_pred_sbm_dyn_decoding(Y, k = k, start = 0, maxit = 200)

# 4. Evaluate clustering accuracy (Adjusted Rand Index)
ari_greedy  <- adjustedRandIndex(fit$cl,  sim_data$U + 1)
ari_viterbi <- adjustedRandIndex(fit$clv, sim_data$U + 1)
cat(sprintf("ARI (Greedy): %.4f | ARI (Viterbi): %.4f\n", ari_greedy, ari_viterbi))

# 5. Model Selection via Dynamic ICL
icl_eval <- icl_dyn_hyb(fit, n = n, TT = TT, k = k)
cat(sprintf("ICL Score: %.2f\n", icl_eval$ICL))

```

### 2. Running the Full Monte Carlo Simulation

To execute the Monte Carlo experiment with perturbed SVD multi-start:

```
source("Simulation.R")

```

Key simulation parameters inside `Simulation.R`:

* `n`: Number of network nodes (default: `20`);

* `TT`: Number of discrete time steps (default: `6`);

* `k`: Number of latent classes (default: `3`);

* `B`: Number of Monte Carlo replicates (default: `100`);

* `nrep`: Number of perturbed spectral restarts per replication (default: `15`);

* `persist`: Markov transition persistence (`"high"` vs. `"low"`).

## 📊 Evaluation Metrics

* **Adjusted Rand Index (ARI)**: Quantifies agreement between recovered community assignments and true ground truth, invariant to label switching.

* **Optimal Label Permutation (`best_perm.R`)**: Identifies the bijective mapping between latent classes to calculate exact classification confusion matrices.

* **Integrated Completed Likelihood (ICL)**: Tailored BIC-like criterion penalizing both static network connectivity parameters and Markov chain transition matrices.
