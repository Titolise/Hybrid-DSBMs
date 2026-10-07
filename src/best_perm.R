require(gtools)

best_perm <- function(Utrue, Uprop) {
  k <- max(Utrue)
  Tab0 <- matrix(0, k, k)
  for (u0 in 1:k) {
    for (u1 in 1:k) {
      Tab0[u0, u1] <- sum(Utrue == u0 & Uprop == u1)
    }
  }
  P <- permutations(k, k, 1:k)
  max_agr <- sum(diag(Tab0))
  perm1 <- 1:k
  Tab1 <- Tab0
  for (j in 1:nrow(P)) {
    Tab <- Tab0[, P[j, ]]
    if (sum(diag(Tab)) > max_agr) {
      perm1 <- P[j, ]
      max_agr <- sum(diag(Tab))
      Tab1 <- Tab
    }
  }
  U1 <- Uprop
  for (j in 1:k) {
    U1[Uprop == perm1[j]] <- j
  }
  out <- list(Tab0 = Tab0, Tab1 = Tab1, perm1 = perm1, U1 = U1)
}
