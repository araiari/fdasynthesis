#' Weighted Karcher Mean of Multivariate Functional Data
#'
#' Calculates (weighted) Karcher mean of a collection of multivariate functional
#' data using the elastic square-root velocity (srvf) framework.
#'
#' @param beta An array of sizes \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#' @param wts Vector of length \eqn{N} with weights to be assigned to each
#'  curve. Use `NULL` for unweighted mean. Defaults to `NULL`.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param maxit Maximum number of iterations.
#' @param ncores XXX
#' @param use_verbose XXX
#'
#' @return Returns a list containing \item{qmu}{mean srvf}
#' \item{betamean}{(Weighted) mean curve}
#' \item{qmu}{(Weighted) mean srvf}
#' \item{beta}{Array of original curves}
#' \item{q}{Array of original srvfs}
#' \item{betan}{Array of aligned curves}
#' \item{qn}{Array of aligned srvfs}
#' \item{wts}{Weights used in the computation}
#' \item{E}{Energy, as the \eqn{L^2}-norm of the mean at each iteration}
#' \item{qun}{Cost function, as the sum of distances to the mean at each iteration}
#' \item{type}{String indicating whether mean or weighted mean is returned}
#'
#' @keywords srvf alignment
#'
#' @importFrom foreach %dopar%
multivariate_weighted_karcher_mean <- function (
    beta,
    wts = NULL,
    lambda = 0.0,
    maxit = 20,
    ncores = 1L,
    use_verbose = TRUE
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

  # Compute number of cores to use
  navail <- max(parallel::detectCores() - 1, 1)

  if (ncores > navail) {
    if (use_verbose)
      cli::cli_alert_warning(
        "The number of requested cores ({ncores}) is larger than the number of
        available cores ({navail}). Using the maximum number of available cores..."
      )
    ncores <- navail
  }

  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else
    foreach::registerDoSEQ()

  # Initialize the mean as the medoid
  v_w = array(sapply(1:N, function(n) q[, , n] * wts[n]), dim = dim(q))
  qmu = rowSums(v_w, dims = 2) / sum(wts)

  n <- NULL
  dists <- foreach::foreach(
    n = 1:N, .combine = "c", .packages = 'fdasrvf'
    ) %dopar% {
    find_rotation_seed_unique(
      qmu, q[ , , n],
      mode = "O",
      alignment = TRUE,
      rotation = FALSE,
      scale = FALSE
    )$d
  }

  medoid_idx <- which.min(dists)
  qmu <- q[, , medoid_idx]
  bmu = beta[, , medoid_idx]

  # Other initializations
  delta = 0.5
  tolv = 1e-04
  told = 5 * 0.001
  itr = 1
  sumd = rep(0, maxit + 1)
  sumd[1] = Inf
  betan = array(0, dim = c(L, M, N))
  qn = array(0, dim = c(L, M, N))
  normbar = rep(0, maxit + 1)

  while (itr < maxit) {

    # align the curves
    alignment_step <- foreach::foreach(
      n = 1:N,
      .combine = cbind,
      .packages = "fdasrvf") %dopar% {
        out <- find_rotation_seed_unique(
          q1 = qmu,
          q2 = q[, , n],
          mode = "O",
          alignment = TRUE,
          rotation = FALSE,
          scale = FALSE,
          lambda = lambda
        )
        list(d = out$d, q2n = out$q2best, gamn = out$gambest)
      }

    d <- unlist(alignment_step[1, ])
    dim(d) <- N
    sumd[itr + 1] <- sum(d^2)

    qn <- unlist(alignment_step[2, ])
    dim(qn) <- c(L, M, N)

    gam <- unlist(alignment_step[3, ])
    dim(gam) <- c(M, N)

    # Mean computation
    v_w = array(sapply(1:N, function(n) qn[, , n] * wts[n]), dim = dim(qn))
    vbar = rowSums(v_w, dims = 2) / sum(wts)
    bbar = fdasrvf::q_to_curve(vbar, scale = 1)

    # # Mean computation alternativa ed equivalente 20240422
    # v_w = array(sapply(1:N, function(n) betan[, , n] * wts[n]), dim = dim(betan))
    # bbar = rowSums(v_w, dims = 2) / sum(wts)
    # vbar = fdasrvf::curve_to_q(bbar, scale = FALSE)$q

    normbar[itr] = sqrt(innerprod_q2(vbar, vbar))


    # Se le distanze dalla media sono aumentate, esci
    if ((sumd[itr] - sumd[itr + 1]) < 0) { break }
    # Se le distanze dalla media sono diminuite a sufficienza e la norma della media non e troppo grande, esci
    if ((normbar[itr] <= tolv) || (abs(sumd[itr + 1] - sumd[itr]) <= told)) { break }

    qmu = vbar
    betamean = bbar

    itr = itr + 1
  }


  # Normalization step
  gam = t(gam)
  gamI = sqrt_weighted_mean_inverse(t(gam))
  betamean = group_action_by_gamma_coord(betamean, gamI) ##
  qmu = fdasrvf::curve_to_q(betamean, scale = FALSE)$q


  ifelse(is_weighted, type <- "Karcher Weighted Median", type <- "Karcher Mean")

  list(
    betamean = betamean,
    qmu = qmu,
    beta = beta,
    q = q,
    betan = betan, ## NULL al momento
    qn = qn,
    wts = wts,
    type = type,
    E = normbar[1:itr],
    qun = sumd[1:itr]
  )

}
