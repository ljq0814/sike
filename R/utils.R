#' Projection Matrix
#'
#' Computes the orthogonal projection matrix onto the column space of a
#' given matrix, defined as \eqn{P_\beta = \beta(\beta^\top \beta)^{-1} \beta^\top}.
#'
#' @param beta A numeric matrix with dimension \code{p x K} and full column rank.
#'
#' @return A \code{p x p} symmetric idempotent projection matrix.
#'
#' @noRd
proj_mat <- function(beta){
  return(beta%*%solve(t(beta)%*%beta)%*%t(beta))
}

#' Trace Correlation between Two Subspaces
#'
#' Computes the trace correlation between two estimated subspaces as a measure
#' of their agreement. A value of 1 indicates perfect alignment.
#'
#' @param beta1 A \code{p x K} matrix representing the estimated subspace.
#' @param beta2 A \code{p x K} matrix representing the true subspace.
#'
#' @return A numeric scalar in \eqn{[0, 1]} giving the trace correlation
#'   between the two subspaces, defined as
#'   \deqn{\frac{1}{K} \text{tr}(P_{\beta_1} P_{\beta_2})}
#'   where \eqn{P_{\beta_1}} and \eqn{P_{\beta_2}} are the projection matrices
#'   onto the column spaces of \code{beta1} and \code{beta2} respectively.
#'
#' @export
trace_corr <- function(beta1, beta2){
  p <- nrow(beta2); K <- ncol(beta2)
  return(sum(diag(proj_mat(beta1)%*%proj_mat(beta2)))/K)
}

