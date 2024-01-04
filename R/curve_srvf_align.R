#' Karcher (weighted) Mean of multidimensional functions
#'
#' Calculates Karcher mean or median of a collection of curves using the elastic
#' square-root velocity (srvf) framework.
#'
#' @param beta Array of sizes \eqn{M \times L \times N} describing \eqn{N}
#' curves of dimension \eqn{L} evaluated on \eqn{M} points
#' @param maxit Maximum number of iterations
#' @param wts Vector of weights for computing the weighted Karcher
#' mean (default = `NULL`)
#'
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{beta}{centered data}
#' \item{betamean}{mean or median curve}
#' \item{type}{string indicating whether mean or median is returned}
#' \item{v}{shooting vectors}
#' \item{q}{array of srvfs}
#' \item{gam}{array of warping functions}
#' \item{cent}{centers of original curves}
#' \item{len}{length of curves}
#' \item{len_q}{length of srvfs}
#' \item{mean_scale}{mean length}
#' \item{mean_scale_q}{mean length srvf}
#' \item{E}{energy}
#' \item{qun}{cost function}
weighted_karcher_mean <- function (beta, lambda = 0.0, maxit = 20, wts = NULL)
{

  if (is.null(wts))
    wts = rep(1, N)

  tmp = dim(beta)
  L = tmp[1]
  M = tmp[2]
  N = tmp[3]
  q = array(0, dim = tmp)
  len = rep(0, N)
  len_q = rep(0, N)
  cent = matrix(0, nrow = L, ncol = N)
  for (ii in 1:N) {
    beta1 = beta[ , , ii]
    centroid1 = calculate_centroid(beta = beta1)
    cent[ , ii] = -1 * centroid1
    dim(centroid1) = c(length(centroid1), 1)
    beta1 = beta1 - repmat(centroid1, 1, M)
    beta[ , , ii] = beta1
    out = fdasrvf::curve_to_q(beta = beta1)
    q[, , ii] = out$q
    len[ii] = out$len
    len_q[ii] = out$lenq
  }

  mu = q[, , 1] # initialize centroid as the first observation
  bmu = beta[, , 1]
  delta = 0.5
  tolv = 1e-04
  told = 5 * 0.001
  itr = 1
  sumd = rep(0, maxit + 1)
  sumd[1] = Inf
  v = array(0, c(L, M, N))
  normvbar = rep(0, maxit + 1)


  cat("\nInitializing...\n")
  gam = matrix(0, M, N)
  for (k in 1:N) {
    out = find_rotation_seed_unqiue(q1 = mu, q2 = q[, , k], mode = "O",
                                    rotated = FALSE, lambda = lambda)
    gam[, k] = out$gambest
  }

  gam = t(gam)
  gamI = fdasrvf::SqrtMeanInverse(gam = t(gam))
  bmu = group_action_by_gamma_coord(f = bmu, gamma = gamI)
  mu = fdasrvf::curve_to_q(beta = bmu)$q # sulla sfera
  mu[is.nan(mu)] = 0

  while (itr < maxit) {
    cat(sprintf("Iteration: %d\n", itr))
    mu = mu / sqrt(innerprod_q2(mu, mu)) # sulla sfera

    for (i in 1:N) {
      q1 = q[ , , i]

      out = find_rotation_seed_unqiue(q1 = mu, q2 = q1, mode = "O",
                                      rotated = FALSE, lambda)
      qn_t = out$q2best / sqrt(innerprod_q2(out$q2best, out$q2best)) #in utils_curve

      q1dotq2 = innerprod_q2(mu, qn_t) #in utils_curve

      if (q1dotq2 > 1){
        q1dotq2 = 1
      }
      if (q1dotq2 < -1){
        q1dotq2 = -1
      }

      dist = acos(q1dotq2)

      u = qn_t - q1dotq2 * q1 # shooting vector di q1 rispetto a qn_t
      normu = sqrt(innerprod_q2(u, u))
      if (normu > 1e-4){
        w = u * acos(q1dotq2)/normu # perché non c'è il seno a denominatore?
      } else {
        w = matrix(0, L, M)
      }

      v[, , i] = w

      sumd[itr + 1] = sumd[itr + 1] + dist^2
    }

    # if(ms == "median"){#run for median only
    #   sumv = rowSums(v_d, dims = 2)
    #   sum_dinv = sum(1/d_i)
    #   vbar = sumv/sum_dinv
    # }
    # else{ #run for mean only
    # sumv = rowSums(v, dims = 2) ##* reso weighted
    # vbar = sumv/N ##*
    v_w = array(sapply(1:2, function(n) v[, , n] * wts[n]), dim = dim(v))
    vbar = rowSums(v_w, dims = 2) / sum(wts)
    # }

    normvbar[itr] = sqrt(innerprod_q2(vbar, vbar))
    normv = normvbar[itr]
    if ((sumd[itr]-sumd[itr+1]) < 0){
      break
    } else if ((normv > tolv) && (abs(sumd[itr + 1] - sumd[itr]) > told)) {
      mu = cos(delta * normvbar[itr]) * mu +
        sin(delta *  normvbar[itr]) * vbar/normvbar[itr] # dallo shooting vector alla media
      x = fdasrvf::q_to_curve(q = mu)
      a = -1 * calculate_centroid(x)
      dim(a) = c(length(a), 1)
      betamean = x + repmat(a, 1, M) #non ci andrebbe un -??
    }
    else {
      break
    }
    itr = itr + 1
  }

  if (scale){
    mean_scale = prod(len)^(1/length(len))
    mean_scale_q = prod(len_q)^(1/length(len))
    betamean = mean_scale*betamean
  }

  ifelse(weighted, type <- "Weighted Karcher mean", type <- "Karcher mean")
  return(list(beta = beta, mu = mu, type = type, betamean = betamean,
              v = v, q = q, E = normvbar[1:itr], cent = cent, len = len,
              len_q = len_q, qun = sumd[1:itr], mean_scale = mean_scale,
              mean_scale_q = mean_scale_q))
}
