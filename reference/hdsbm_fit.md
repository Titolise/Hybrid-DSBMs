# Hybrid DSBM Estimation with Perturbed SVD Multi-Start

Estimates a Dynamic Stochastic Block Model (DSBM) by combining the
Hybrid CEM algorithm with a multi-start strategy based on the spectral
decomposition (SVD) of the time-averaged adjacency matrix and stochastic
perturbations.

## Usage

``` r
hdsbm_fit(
  Y,
  k,
  nrep = 10,
  perturb_range = c(0.05, 0.35),
  maxit = 200,
  tol = 1e-10,
  verbose = TRUE,
  seed = NULL
)
```

## Arguments

- Y:

  Binary 3D array of dimensions \\n \times n \times T\\ (the temporal
  adjacency matrix).

- k:

  Integer, number of latent communities/blocks.

- nrep:

  Integer, number of random restarts perturbed via SVD (default: 10). If
  0, only the initial fit is performed.

- perturb_range:

  Numeric vector of length 2 indicating the lower and upper limits of
  the random perturbation rate (default: `c(0.05, 0.35)`).

- maxit:

  Integer, maximum number of CEM iterations per restart (default: 200).

- tol:

  Numeric, convergence tolerance for the log-likelihood (default:
  1e-10).

- verbose:

  Logical; if `TRUE`, prints restart progress to the console (default:
  `TRUE`).

- seed:

  Optional integer to ensure reproducibility.

## Value

An object of class `hdsbm` containing:

- `best_fit`: the estimate with the highest log-likelihood (parameters
  \\\pi_0, \Pi, \Psi\\, classifications `cl` and `clv`).

- `lk_trace`: vector containing the log-likelihoods obtained from each
  restart. \#'

- `best_restart`: the index of the restart that achieved the optimal
  value (0 = baseline).

- `icl`: ICL score calculated on the optimal fit.

- `params`: model dimensional parameters (`n`, `k`, `TT`).
