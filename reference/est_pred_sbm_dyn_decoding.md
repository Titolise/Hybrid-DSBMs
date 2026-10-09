# Estimation of the Hybrid DSBM model via the CEM algorithm

Estimation of the Hybrid DSBM model via the CEM algorithm

## Usage

``` r
est_pred_sbm_dyn_decoding(
  Y,
  k,
  start = 0,
  tol = 10^-10,
  maxit = 100,
  Tau = NULL
)
```

## Arguments

- Y:

  3D binary array (n x n x TT) representing the temporal adjacency
  matrix.

- k:

  Number of latent blocks/classes.

- start:

  Initialization strategy (0 = K-means, 1 = random, 2 = custom Tau).

- tol:

  Convergence tolerance for the log-likelihood.

- maxit:

  Maximum number of CEM iterations.

- Tau:

  Initial membership array (required if start = 2).

## Value

List containing parameter estimates, likelihood, and classifications (cl
and clv).
