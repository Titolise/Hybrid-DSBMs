#' Synthetic dataset for Dynamic SBM models
#'
#' A sample dataset containing a dynamic network simulated under
#' the Hybrid DSBM framework.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{Y}{Three-dimensional binary array (25 x 25 x 4).}
#'   \item{U_true}{25 x 4 matrix of true community memberships.}
#'   \item{params}{Dimensional parameters (n, k, TT).}
#' }
#' @source Data simulated via \code{draw_sn_dyn}.
#' @examples
#' data(toy_hdsbm)
#' dim(toy_hdsbm$Y)
"toy_hdsbm"