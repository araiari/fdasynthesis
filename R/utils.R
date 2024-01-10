#' Estimate w_a
#'
#' @description
#' Estimate the optimal value of `w_a` for convex combination of two distances.
#'
#' @details
#' This functions finds the value of `w_a` to combine two distance matrices
#' `D1` and `D2`. The optimal value is set as the one maximizing the
#' compromise between hierarchical clustering output with `D1` and with `D2`.
#' It follows from the algorithm `hclustcompro` in the library `SPARTAAS`.
#'
#' @param D1 matrix of dimension \eqn{N \times N} : first distance matrix
#' @param D2 matrix of dimension \eqn{N \times N} : second distance matrix
#' @param n_iterations integer in \eqn{\{1,\dots,N-1\}} : number of iteration for
#' the choice of alpha. Default is set to 5.
#' @return `w_a` : scalar, with \eqn{w_a \in [0,1]} : estimated optimal parameter.
estimate_optimal_w_a = function(
    D1, D2, n_iterations = 5
){
  out = SPARTAAS::hclustcompro_select_alpha(
    D1 = D1,
    D2 = D2,
    acc = 2,
    resampling = T,
    method = "complete",
    iter = n_iterations
  )
  out$alpha
}


#' Compute the weights of each of the neighbors for weighted averaging
#'
#' @description
#' This functions computes the weights to be assigned to each neighbors:
#' it looks at the distances of the k closest elements, sample weights
#' from a Dirichlet of parameters proportional to such distances, return
#' indices of the neighbors and weights.
#'
#' @param D_i vector of length \eqn{N} : distances between one item and all the others
#' @param K scalar : number of neighbors to be considered.
#' @param alpha_0  scalar, with \eqn{\alpha_0 >0} : concentation parameter of the
#' Dirichlet distribution from which the weights of neighbors are sampled.
#' @param beta_0 scalar, with \eqn{\beta_0 > 0} : rate parameter
#' weighting the inverse distances in the weight computation
#' @return: list of two elements:
#' - `index_k` vector of length `K` : indices of the `K` closest neighbors
#' - `p_k`  vector of length `K` : computed weights, one for each neighbors
compute_weights = function(
    D_i, K, alpha_0 = 1, beta_0 = 1
){
  # Finding the k closest curves within the cluster
  index_k = order(D_i)[1+1:K]

  # From the distances, compute the weights p_k
  d_k = D_i[index_k]

  inv_d_k = exp(-beta_0 * d_k)
  alpha_k = alpha_0*inv_d_k / sum(inv_d_k)

  p_k = c(rdirichlet(n=1, alpha=alpha_k))

  list(index_k=index_k, p_k=p_k)
}



#' Random sample from a Dirichlet distribution
#'
#' @param n integer : number of items to be sampled
#' @param alpha scalar, with \eqn{\alpha > 0} : concentration parameter of the
#' distribution
#' @param method method for sampling. Either `"beta"` or `"gamma"`.
#'
#' @return A numeric matrix of shape \eqn{n \times \mathrm{length}(\alpha)}
#'   storing an $n$-sample drawn from a Dirichlet distribution with parameters
#'   \eqn{\alpha}.
rdirichlet <- function(n, alpha, method = c("beta", "gamma")) {
  method <- rlang::arg_match(method)
  switch(
    method,
    beta = rdirichlet_via_beta(n = n, alpha = alpha),
    gamma = rdirichlet_via_gamma(n = n, alpha = alpha)
  )
}



#' Random sample from a Dirichlet distribution through gamma
#'
#' @param n integer : number of items to be sampled
#' @param alpha scalar, with \eqn{\alpha > 0} : concentration parameter of the
#' distribution
#'
#' @return A numeric matrix of shape \eqn{n \times \mathrm{length}(\alpha)}
#'   storing an $n$-sample drawn from a Dirichlet distribution with parameters
#'   \eqn{\alpha}.
rdirichlet_via_gamma <- function(n, alpha) {
  y <- lapply(alpha, function(.alpha) stats::rgamma(n = n, shape = .alpha, rate = 1))
  y <- do.call(cbind, y)
  s <- rowSums(y)
  n_nan <- sum(s < .Machine$double.eps)
  while (n_nan > 0) {
    yl <- lapply(alpha, function(.alpha) stats::rgamma(n = n_nan, shape = .alpha, rate = 1))
    yl <- do.call(cbind, yl)
    y[s < .Machine$double.eps, ] <- yl
    s <- rowSums(y)
    n_nan <- sum(s < .Machine$double.eps)
  }
  s <- matrix(s, nrow = n, ncol = length(alpha))
  y / s
}


#' Random sample from a Dirichlet distribution through beta
#'
#' @param n integer : number of items to be sampled
#' @param alpha scalar, with \eqn{\alpha > 0} : concentration parameter of the
#' distribution
#'
#' @return A numeric matrix of shape \eqn{n \times \mathrm{length}(\alpha)}
#'   storing an $n$-sample drawn from a Dirichlet distribution with parameters
#'   \eqn{\alpha}.
rdirichlet_via_beta <- function(n, alpha) {
  D <- length(alpha)
  if (D == 1) return(matrix(rep(1, n), nrow = n))
  x <- list()
  x[[1]] <- stats::rbeta(n = n, shape1 = alpha[1], shape2 = sum(alpha[-1]))
  if (D > 2) {
    for (j in 2:(D - 1)) {
      phi <- stats::rbeta(n = n, shape1 = alpha[j], shape2 = sum(alpha[(j + 1):D]))
      x_prev <- x[1:(j - 1)]
      x_prev <- do.call(cbind, x_prev)
      x_prev <- rowSums(x_prev)
      x[[j]] <- (1 - x_prev) * phi
    }
  }
  x_prev <- x[1:(D - 1)]
  x_prev <- do.call(cbind, x_prev)
  x_prev <- rowSums(x_prev)
  x[[D]] <- 1 - x_prev
  do.call(cbind, x)
}



#' Integration
#'
#' @param x vector
#' @param y function
#'
#' @return scalar
integrate = function (
    x, y
) {
  integrand = stats::splinefun(x, y)
  stats::integrate(
    f = integrand, lower = min(x), upper = max(x)
    )$value
}



#' Cumulative integration
#'
#' @param x vector
#' @param y function
#'
#' @return vector
cumintegrate = function (
    x, y
) {
  inegrand = stats::splinefun(x, y)
  purrr::map_dbl(x, \(.x) {
    stats::integrate(
      f = inegrand,
      lower = min(x),
      upper = .x
    )$value
  })
}

linear_index <- function(n) {
  res <- tidyr::expand_grid(i = 1:n, j = 1:n)
  res <- subset(res, res$j > res$i)
  res
}
