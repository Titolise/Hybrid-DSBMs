# Synthetic dataset for Dynamic SBM models

A sample dataset containing a dynamic network simulated under the Hybrid
DSBM framework.

## Usage

``` r
toy_hdsbm
```

## Format

A list with 3 elements:

- Y:

  Three-dimensional binary array (25 x 25 x 4).

- U_true:

  25 x 4 matrix of true community memberships.

- params:

  Dimensional parameters (n, k, TT).

## Source

Data simulated via `draw_sn_dyn`.

## Examples

``` r
data(toy_hdsbm)
dim(toy_hdsbm$Y)
#> [1] 25 25  4
```
