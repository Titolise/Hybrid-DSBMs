#' Simulation of Dynamic Networks under the DSBM Model
#'
#' Generates a temporal adjacency tensor and latent membership trajectories
#' for observed nodes across discrete time windows.
#'
#' @param n Integer, number of nodes.
#' @param k Integer, number of latent classes/communities.
#' @param TT Integer, number of observed time windows.
#' @param piv Numeric vector of size \eqn{k}, prior probabilities for the initial time point.
#' @param Pi Markov transition matrix of size \eqn{k \times k}.
#' @param Psi Matrix of block connection probabilities (\eqn{k \times k}).
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{U}: Matrix of size \eqn{n \times TT} containing the true latent classes (0-indexed).
#'   \item \code{Y}: Binary 3D array of size \eqn{n \times n \times TT}.
#' }
#' @importFrom stats rmultinom runif
#' @export

draw_sn_dyn <- function(n, k, TT, piv, Pi, Psi) {
  U <- matrix(0, n, TT)
  if (k > 1) {
    for (i in 1:n) {
      U[i, 1] <- which(stats::rmultinom(1, 1, piv) == 1) - 1
      for (t in 2:TT) {
        U[i, t] <- which(stats::rmultinom(1, 1, Pi[U[i, t - 1] + 1, ]) == 1) - 1
      }
    }
  }
  
  Y <- array(0, c(n, n, TT))
  for (t in 1:TT) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (k == 1) {
          Y[i, j, t] <- Y[j, i, t] <- 1 * (stats::runif(1) < Psi)
        } else {
          Y[i, j, t] <- Y[j, i, t] <- 1 * (stats::runif(1) < Psi[U[i, t] + 1, U[j, t] + 1])
        }
      }
    }
  }
  
  list(U = U, Y = Y)
}