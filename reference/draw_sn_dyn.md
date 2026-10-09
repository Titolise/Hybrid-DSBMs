# Simulation of Dynamic Networks under the DSBM Model

Generates a temporal adjacency tensor and latent membership trajectories
for observed nodes across discrete time windows.

## Usage

``` r
draw_sn_dyn(n, k, TT, piv, Pi, Psi)
```

## Arguments

- n:

  Integer, number of nodes.

- k:

  Integer, number of latent classes/communities.

- TT:

  Integer, number of observed time windows.

- piv:

  Numeric vector of size \\k\\, prior probabilities for the initial time
  point.

- Pi:

  Markov transition matrix of size \\k \times k\\.

- Psi:

  Matrix of block connection probabilities (\\k \times k\\).

## Value

A list with two components:

- `U`: Matrix of size \\n \times TT\\ containing the true latent classes
  (0-indexed).

- `Y`: Binary 3D array of size \\n \times n \times TT\\.
