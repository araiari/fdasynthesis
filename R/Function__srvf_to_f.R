#' Transformation from SRSF Space
#'
#' This function transforms SRVFs back to the original functional space.
#'
#' @param q Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the SRSFs that need to be transformed.
#'   
#'   - If a vector, it must be of shape \eqn{M} and it is interpreted as a
#'   single \eqn{1}-dimensional curve observed on a grid of size \eqn{M}.
#'   - If a matrix and `multidimensional == FALSE`, it must be of shape
#'   \eqn{M \times N}. In this case, it is interpreted as a sample of \eqn{N}
#'   curves observed on a grid of size \eqn{M}, unless \eqn{M = 1} in which case
#'   it is interpreted as a single \eqn{1}-dimensional curve observed on a grid
#'   of size \eqn{M}.
#'   - If a matrix and `multidimensional == TRUE`,it must be of shape
#'   \eqn{L \times M} and it is interpreted as a single \eqn{L}-dimensional
#'   curve observed on a grid of size \eqn{M}.
#'   - If a 3D array, it must be of shape \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   SRSFs are evaluated.
#' @param f0 Either a numeric value or a numeric vector of or a numeric matrix
#'   specifying the initial value of the curves in the original functional
#'   space. It must be:
#'   
#'   - a value if `q` represents a single \eqn{1}-dimensional SRSF.
#'   - a vector of length \eqn{L} if `q` represents a single
#'   \eqn{L}-dimensional SRSF.
#'   - a vector of length \eqn{N} if `q` represents a sample of \eqn{N}
#'   \eqn{1}-dimensional SRSFs.
#'   - a matrix of shape \eqn{L \times M} if `q` represents a sample of \eqn{N}
#'   \eqn{L}-dimensional  SRSFs.
#' @param multidimensional A boolean specifying if the curves are
#'   multi-dimensional. This is useful when `q` is provided as a matrix to
#'   determine whether it is a single multi-dimensional curve or a collection of
#'   uni-dimensional curves. Defaults to `FALSE`.
#'
#' @return A numeric array of the same shape as the input `q` storing the
#'   transformation of the SRSFs `q` back to the original functional space.

srvf_to_f = function (q, time, f0 = NULL, multidimensional = FALSE) 
{
  dims <- dim(q)
  dims0 <- dim(f0)
  if (is.null(dims)) {
    L <- 1
    stopifnot(length(f0) == L)
    M <- length(q)
    N <- 1
    integrand <- q * abs(q)
    if (is.null(f0)) f0 = 0 ##*
    f <- f0 + pracma::cumtrapz(time, integrand)
  }
  else {
    if (length(dims) == 2) {
      if (multidimensional || dims[1] == 1) {
        L <- dims[1]
        stopifnot(is.null(dims0) && length(f0) == L)
        M <- dims[2]
        N <- 1
        norm_q <- sqrt(colSums(q^2))
        if (is.null(f0)) f0 = rep(0, L) ##*
        f <- lapply(1:L, function(l) {
          integrand <- q[l, ] * norm_q
          f0[l] + pracma::cumtrapz(time, integrand)
        })
        f_temp <- do.call(rbind, f) ################ da cambiare
        f = matrix(nrow=L, ncol=M) ##*
        for (l in 1:L) f[l,] = f_temp[1:M+M*(l-1)] ##*
      }
      else {
        L <- 1
        M <- dims[1]
        N <- dims[2]
        stopifnot(is.null(dims0) && length(f0) == N)
        if (is.null(f0)) f0 = rep(0, N) ##*
        f <- lapply(1:N, function(n) {
          integrand <- q[, n] * abs(q[, n])
          f0[n] + pracma::cumtrapz(time, integrand)
        })
        f <- do.call(cbind, f)
      }
    }
    else {
      stopifnot(!is.null(dims0) && length(dims0) == 2)
      L <- dims[1]
      stopifnot(dims0[1] == L)
      M <- dims[2]
      N <- dims[3]
      stopifnot(dims0[2] == N)
      norm_q <- lapply(1:N, function(n) {
        sqrt(colSums(q[, , n]^2))
      })
      if (is.null(f0)) f0 = matrix(0, nrow=L, ncol=N) ##*
      f <- lapply(1:N, function(n) {
        res <- lapply(1:L, function(l) {
          integrand <- q[l, , n] * norm_q[[n]][l ] ##*
          f0[l, n] + pracma::cumtrapz(time, integrand)
        })
        do.call(rbind, res)
      })
      f_temp <- do.call(cbind, f)
      f = array(dim = c(L,M,N)) ##*
      for (l in 1:L) f[l,,] = f_temp[1:M+M*(l-1),] ##*
    }
  }
  f
}
