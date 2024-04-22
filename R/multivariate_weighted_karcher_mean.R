#' Weighted Karcher Mean of Multivariate Functional Data
#'
#' Calculates (weighted) Karcher mean of a collection of multivariate functional
#' data using the elastic square-root velocity (srvf) framework.
#'
#' @param beta An array of sizes \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param maxit Maximum number of iterations.
#' @param wts Vector of length \eqn{N} with weights to be assigned to each
#'  curve. Use `NULL` for unweighted mean. Defaults to `NULL`.
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{betamean}{(Weighted) mean curve}
#' \item{mu}{(Weighted) mean srvf}
#' \item{beta}{Array of original curves}
#' \item{q}{Array of original srvfs}
#' \item{betan}{Array of aligned curves}
#' \item{qn}{Array of aligned srvfs}
#' \item{wts}{Weights used in the computation}
#' \item{E}{Energy, as the \eqn{L^2}-norm of the mean at each iteration}
#' \item{qun}{Cost function, as the sum of distances to the mean at each iteration}
#' \item{type}{String indicating whether mean or weighted mean is returned}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'    analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine
#'    Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- multivariate_weighted_karcher_mean(fdasrvf::beta[, , 1, 1:2], maxit = 2)
#' # note: use more functions, small for speed
multivariate_weighted_karcher_mean <- function (
    beta,
    lambda = 0.0,
    maxit = 20,
    wts = NULL
) {

  tmp = dim(beta)
  L = tmp[1]
  M = tmp[2]
  N = tmp[3]

  is_weighted = !(is.null(wts))
  if (!is_weighted)
    wts = rep(1, N)

  time1 <- seq(0, 1, length.out = M)
  q = array(0, dim = c(L, M, N))
  for (ii in 1:N) {
    out = fdasrvf::curve_to_q(beta[, , ii], scale = FALSE)
    q[, , ii] = out$q
  }

  mu = q[, , 1]
  bmu = beta[, , 1]
  delta = 0.5
  tolv = 1e-04
  told = 5 * 0.001
  itr = 1
  sumd = rep(0, maxit + 1)
  sumd[1] = Inf
  betan = array(0, dim = c(L, M, N))
  qn = array(0, dim = c(L, M, N))
  normbar = rep(0, maxit + 1)


  cat("\nInitializing...\n")
  gam = matrix(0, nrow = M, ncol = N)
  for (k in 1:N) {
    out = find_rotation_seed_unqiue(mu, q[, , k], lambda) ##
    gam[, k] = out$gambest
  }

  gamI = SqrtWeightedMeanInverse(gam) ##
  bmu = group_action_by_gamma_coord(bmu, gamI) ##
  mu = fdasrvf::curve_to_q(bmu, FALSE)$q
  mu[is.nan(mu)] = 0

  while (itr < maxit) {
    cat(sprintf("Iteration: %d\n", itr))
    mu = mu

    for (i in 1:N) {
      q1 = q[, , i]

      out = find_rotation_seed_unqiue(mu, q1, lambda) ##
      dist = sqrt( sum((mu - out$q2best)^2) / M)

      qn[, , i] = out$q2best
      betan[, , i] = fdasrvf::q_to_curve(out$q2best, scale = 1)

      gam[,i] = out$gambest

      sumd[itr + 1] = sumd[itr + 1] + dist^2
    }
    cat("\nQui1\n")

    # Mean computation
    v_w = array(sapply(1:N, function(n) qn[, , n] * wts[n]), dim = dim(qn))
    cat("\nQui2\n")
    vbar = rowSums(v_w, dims = 2) / sum(wts)
    cat("\nQui3\n")
    bbar = fdasrvf::q_to_curve(vbar, scale = 1)

    cat("\nQui4\n")
    normbar[itr] = sqrt(innerprod_q2(vbar, vbar))

    # Se le distanze dalla media sono aumentate, esci
    if ((sumd[itr] - sumd[itr + 1]) < 0) { break }
    # Se le distanze dalla media sono diminuite a sufficienza e la norma della media non e troppo grande, esci
    if ((normbar[itr] <= tolv) || (abs(sumd[itr + 1] - sumd[itr]) <= told)) { break }

    mu = vbar
    betamean = bbar

    itr = itr + 1
  }


  # Normalization step
  gam = t(gam)
  gamI = SqrtWeightedMeanInverse(t(gam)) ##
  betamean = group_action_by_gamma_coord(betamean, gamI) ##
  mu = fdasrvf::curve_to_q(betamean, scale = FALSE)$q


  ifelse(is_weighted, type <- "Karcher Weighted Median", type <- "Karcher Mean")
  return(list(betamean = betamean, mu = mu, beta = beta, q = q, betan = betan,
              qn = qn,  wts = wts, type = type, E = normbar[1:itr],
              qun = sumd[1:itr]))
}
