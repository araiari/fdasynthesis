#' SRVF transform of warping functions
#'
#' This function calculates the srvf of warping functions with corresponding
#' shooting vectors and finds the inverse of (weighted) mean.
#'
#' It is based on the `SqrtMeanInverse` funcion of the `fdasrvf` package, and
#' enlarges it potentiality by computing the mean both non-weighted and weighted
#'
#' @param gam matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with
#' \eqn{N} samples
#' @param wts vector of weights
#' @return `gamI` inverse of (weighted) Karcher mean warping function
#'
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' gamI <- SqrtWeightedMeanInverse(fdasrvf::simu_warp$warping_functions)
SqrtWeightedMeanInverse <- function(
    gam,
    wts = NULL
){

  TT = nrow(gam)
  n = ncol(gam)
  eps = .Machine$double.eps
  time <- seq(0, 1, length.out=TT)

  if(is.null(wts))
    wts = rep(1, n)

  psi = matrix(0, TT, n)
  binsize <- mean(diff(time))
  for (i in 1:n){
    psi[,i] = sqrt(fdasrvf::gradient(gam[, i],binsize))
  }

  # Find Direction for initialization
  mu = rowSums(psi * matrix(wts, nrow=TT, ncol=n, byrow=T)) / sum(wts)

  stp <- .3
  maxiter = 501
  vec = matrix(0,TT,n)
  lvm = rep(0,maxiter)
  iter <- 1

  for (i in 1:n){
    vec[,i] <- inv_exp_map(mu, psi[,i])
  }
  vbar <- rowMeans(vec)
  lvm[iter] <- l2_norm(vbar)

  while (lvm[iter]>0.00000001 & iter<maxiter){
    mu <- exp_map(mu, stp*vbar)
    iter <- iter + 1
    for (i in 1:n){
      vec[,i] <- inv_exp_map(mu, psi[,i])
    }
    vbar <- rowSums(vec * matrix(wts, nrow=TT, ncol=n, byrow=T)) / sum(wts)
    lvm[iter] <- l2_norm(vbar)
  }

  gam_mu = cumintegrate(time, mu*mu)
  gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))
  gamI = fdasrvf::invertGamma(gam_mu)
  return(gamI)
}
