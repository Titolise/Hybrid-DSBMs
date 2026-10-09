# Integrated Completed Likelihood (ICL) for Dynamic SBM

Computes the Integrated Completed Likelihood (ICL) criterion for Dynamic
Stochastic Block Models. Penalizes dyadic connection probabilities,
first-order Markov transition parameters, and initial community
proportions.

## Usage

``` r
icl(object, n = NULL, TT = NULL, k = NULL, detailed = FALSE, ...)

# S3 method for class 'hdsbm'
icl(object, n = NULL, TT = NULL, k = NULL, detailed = FALSE, ...)

# Default S3 method
icl(object, n = NULL, TT = NULL, k = NULL, detailed = FALSE, ...)

ICL(object, ...)

icl_dyn_hyb(out_pred, n = NULL, TT = NULL, k = NULL)
```

## Arguments

- object:

  An object of class `hdsbm` (returned by
  [`hdsbm_fit`](https://titolise.github.io/Hybrid-DSBMs/reference/hdsbm_fit.md)),
  or a list returned by
  [`est_pred_sbm_dyn_decoding`](https://titolise.github.io/Hybrid-DSBMs/reference/est_pred_sbm_dyn_decoding.md).

- n:

  Optional integer, number of nodes. Inferred automatically if `NULL`.

- TT:

  Optional integer, number of time points. Inferred automatically if
  `NULL`.

- k:

  Optional integer, number of latent classes. Inferred automatically if
  `NULL`.

- detailed:

  Logical; if `FALSE` (default), returns a single numeric value (like
  `AIC` or `BIC`). If `TRUE`, returns a list with detailed penalty
  decompositions.

- ...:

  Additional arguments passed to methods.

- out_pred:

  Alias for `object` for backward compatibility.

## Value

If `detailed = FALSE`, a numeric value representing the ICL score. If
`detailed = TRUE`, a list containing:

- `ICL`: The total ICL score.

- `lk`: The complete-data log-likelihood.

- `pen_total`: The sum of all penalties.

- `pen_components`: Vector with individual penalties (`psi`, `rho`,
  `pi`).

## Examples

``` r
data(toy_hdsbm)
model <- hdsbm_fit(toy_hdsbm$Y, k = 3, nrep = 2, maxit = 20, verbose = FALSE)
#> ------------|-------------|-------------|-------------|
#>   iteration |   classes   |      lk     |    lk-lko   |
#> ------------|-------------|-------------|-------------|
#>           0 |           3 |    -998.141 |          NA | 
#>           1 |           3 |    -957.418 |     40.7225 | 
#>           2 |           3 |    -949.153 |     8.26474 | 
#>           3 |           3 |    -948.342 |    0.811689 | 
#>           4 |           3 |    -948.278 |   0.0633762 | 
#>           5 |           3 |    -947.982 |    0.296342 | 
#>           6 |           3 |    -947.833 |    0.148739 | 
#>           7 |           3 |    -947.814 |   0.0194721 | 
#>           8 |           3 |    -947.811 |   0.0031965 | 
#>           9 |           3 |     -947.81 | 0.000621065 | 
#>          10 |           3 |     -947.81 | 0.000147667 | 
#>          11 |           3 |     -947.81 | 4.19949e-05 | 
#>          12 |           3 |     -947.81 | 1.34213e-05 | 
#>          13 |           3 |     -947.81 | 4.55998e-06 | 
#>          14 |           3 |     -947.81 |  1.5941e-06 | 
#>          15 |           3 |     -947.81 | 5.64384e-07 | 
#>          16 |           3 |     -947.81 | 2.00939e-07 | 
#>          17 |           3 |     -947.81 | 7.17207e-08 | 
#>          17 |           3 |     -947.81 | 7.17207e-08 | 
#> ------------|-------------|-------------|-------------|
#> ------------|-------------|-------------|-------------|
#>   iteration |   classes   |      lk     |    lk-lko   |
#> ------------|-------------|-------------|-------------|
#>           0 |           3 |    -1114.28 |          NA | 
#>           1 |           3 |    -995.667 |     118.614 | 
#>           2 |           3 |    -950.136 |     45.5307 | 
#>           3 |           3 |    -949.411 |    0.725454 | 
#>           4 |           3 |    -947.712 |     1.69872 | 
#>           5 |           3 |    -947.422 |    0.289896 | 
#>           6 |           3 |    -947.411 |   0.0108429 | 
#>           7 |           3 |    -947.411 | 0.000651327 | 
#>           8 |           3 |    -947.411 | 4.47358e-05 | 
#>           9 |           3 |    -947.411 | 3.55169e-06 | 
#>          10 |           3 |    -947.411 | 3.32586e-07 | 
#>          11 |           3 |    -947.411 | 3.60264e-08 | 
#>          11 |           3 |    -947.411 | 3.60264e-08 | 
#> ------------|-------------|-------------|-------------|
#> ------------|-------------|-------------|-------------|
#>   iteration |   classes   |      lk     |    lk-lko   |
#> ------------|-------------|-------------|-------------|
#>           0 |           3 |    -1145.27 |          NA | 
#>           1 |           3 |     -1041.4 |     103.861 | 
#>           2 |           3 |     -956.86 |     84.5447 | 
#>           3 |           3 |    -950.427 |     6.43343 | 
#>           4 |           3 |     -950.33 |   0.0965194 | 
#>           5 |           3 |    -950.325 |  0.00539881 | 
#>           6 |           3 |    -950.324 | 0.000699781 | 
#>           7 |           3 |    -950.324 | 0.000120465 | 
#>           8 |           3 |    -950.324 | 2.44019e-05 | 
#>           9 |           3 |    -950.324 | 5.42853e-06 | 
#>          10 |           3 |    -950.324 | 1.27035e-06 | 
#>          11 |           3 |    -950.324 | 3.05376e-07 | 
#>          12 |           3 |    -950.324 |  7.4475e-08 | 
#>          12 |           3 |    -950.324 |  7.4475e-08 | 
#> ------------|-------------|-------------|-------------|

# Standard numeric value (ideal for model comparison)
icl(model)
#> [1] -984.8522

# Detailed penalty breakdown
icl(model, detailed = TRUE)
#> $ICL
#> [1] -984.8522
#> 
#> $lk
#> [1] -947.4106
#> 
#> $pen_total
#> [1] 37.44157
#> 
#> $pen_components
#>       psi       rho        pi 
#> 21.270231 12.952464  3.218876 
#> 
```
