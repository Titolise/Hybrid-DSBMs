## code to prepare `toy_hdsbm` dataset goes here

## Preparation code for toy_hdsbm
set.seed(123)

n <- 25
k <- 3
TT <- 4

piv0 <- rep(1 / k, k)
Pi0 <- matrix(0.1, k, k)
diag(Pi0) <- 0.8
Pi0 <- sweep(Pi0, 1, rowSums(Pi0), "/")

Psi0 <- matrix(0.05, k, k)
diag(Psi0) <- 0.45

# Generation with package internal functions
sim <- draw_sn_dyn(n = n, k = k, TT = TT, piv = piv0, Pi = Pi0, Psi = Psi0)

toy_hdsbm <- list( 
  Y = sim$Y, 
  U_true = sim$U + 1, 
  params = list(n = n, k = k, TT = TT)
)

# Save the compressed object to data/toy_hdsbm.rda
usethis::use_data(toy_hdsbm, overwrite = TRUE, compress = "xz")